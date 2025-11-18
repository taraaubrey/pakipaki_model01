"""
Analyze PEST++ IES phi results to understand arr-h observation group performance
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Set paths
master_dir = r'models/local_run33/pest/master_ies'

# Read phi group data
phi_group_file = os.path.join(master_dir, 'local_run33.phi.group.csv')
phi_group = pd.read_csv(phi_group_file)

print("="*80)
print("PHI GROUP ANALYSIS")
print("="*80)
print(f"\nPhi group dataframe shape: {phi_group.shape}")
print(f"Columns: {phi_group.columns.tolist()}")
print("\nFirst few rows:")
print(phi_group.head(10))

# Get observation group columns (excluding iteration and real_name)
# Filter out non-numeric columns
metadata_cols = ['iteration', 'total_runs', 'obs_realization', 'par_realization', 'real_name', 'total_phi']
obs_groups = [col for col in phi_group.columns if col not in metadata_cols]
print(f"\nObservation groups found: {obs_groups}")

# Convert to numeric (some columns might have mixed types)
for col in obs_groups:
    phi_group[col] = pd.to_numeric(phi_group[col], errors='coerce')

# Calculate mean phi by iteration for each group
phi_by_iter = phi_group.groupby('iteration')[obs_groups].mean()
print("\n" + "="*80)
print("MEAN PHI BY ITERATION FOR EACH OBSERVATION GROUP")
print("="*80)
print(phi_by_iter)

# Calculate phi reduction from iteration 0
print("\n" + "="*80)
print("PHI REDUCTION (%) FROM ITERATION 0")
print("="*80)
phi_reduction = (phi_by_iter.iloc[0] - phi_by_iter) / phi_by_iter.iloc[0] * 100
print(phi_reduction)

# Identify arr-h related groups
arr_h_groups = [col for col in obs_groups if 'arr' in col.lower() or 'obg' in col.lower()]
print(f"\n" + "="*80)
print(f"POTENTIAL ARR-H GROUPS: {arr_h_groups}")
print("="*80)

# Read phi_factors to understand group weights
phi_factors_file = os.path.join(master_dir, 'phi_factors.csv')
phi_factors = pd.read_csv(phi_factors_file, header=None, names=['group', 'weight'])
print("\n" + "="*80)
print("PHI FACTORS (WEIGHTS)")
print("="*80)
print(phi_factors)

# Read obs_data to map observation names to groups (read first 1000 rows)
obs_data_file = os.path.join(master_dir, 'local_run33.obs_data.csv')
print("\nReading observation data (first 10000 rows)...")
obs_data = pd.read_csv(obs_data_file, nrows=10000)
print(f"\nObs data shape: {obs_data.shape}")
print(f"Columns: {obs_data.columns.tolist()}")
print("\nFirst few rows:")
print(obs_data.head(20))

# Find arr-h observations
if 'obgnme' in obs_data.columns and 'oname' in obs_data.columns:
    arr_h_obs = obs_data[obs_data['oname'].str.contains('arr', case=False, na=False)]
    print("\n" + "="*80)
    print(f"ARR-H OBSERVATIONS (sample from first 10k):")
    print("="*80)
    print(arr_h_obs[['oname', 'obgnme', 'obsval']].head(20))

    print("\nUnique observation groups for arr-h observations:")
    print(arr_h_obs['obgnme'].unique())

# Calculate final vs initial phi for each group
print("\n" + "="*80)
print("PHI PERFORMANCE SUMMARY")
print("="*80)
final_iter = phi_by_iter.index[-1]
performance = pd.DataFrame({
    'Initial_Phi': phi_by_iter.iloc[0],
    'Final_Phi': phi_by_iter.iloc[-1],
    'Reduction_%': ((phi_by_iter.iloc[0] - phi_by_iter.iloc[-1]) / phi_by_iter.iloc[0] * 100),
    'Ratio_Final/Initial': phi_by_iter.iloc[-1] / phi_by_iter.iloc[0]
})
performance = performance.sort_values('Reduction_%', ascending=False)
print(performance)

# Plot phi evolution
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

# Plot 1: All groups over iterations
ax = axes[0]
for col in obs_groups:
    ax.plot(phi_by_iter.index, phi_by_iter[col], marker='o', label=col)
ax.set_xlabel('Iteration')
ax.set_ylabel('Mean Phi')
ax.set_title('Phi Evolution by Observation Group')
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
ax.grid(True, alpha=0.3)

# Plot 2: Phi reduction percentage
ax = axes[1]
phi_reduction.T.plot(ax=ax, kind='bar', width=0.8)
ax.set_xlabel('Observation Group')
ax.set_ylabel('Phi Reduction (%)')
ax.set_title('Phi Reduction from Initial (Iteration 0)')
ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax.grid(True, alpha=0.3, axis='y')
ax.legend(title='Iteration', bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

plt.tight_layout()
plt.savefig('models/local_run33/phi_analysis.png', dpi=150, bbox_inches='tight')
print("\n" + "="*80)
print("Figure saved to: models/local_run33/phi_analysis.png")
print("="*80)

# Additional analysis: contribution to total phi
print("\n" + "="*80)
print("CONTRIBUTION TO TOTAL PHI (% of total)")
print("="*80)
phi_contribution = phi_by_iter.div(phi_by_iter.sum(axis=1), axis=0) * 100
print(phi_contribution)

print("\nContribution at final iteration:")
print(phi_contribution.iloc[-1].sort_values(ascending=False))
