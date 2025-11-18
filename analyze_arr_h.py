"""
Detailed analysis of arr-h observation groups
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

# Define arr-h groups based on the analysis
arr_h_groups = {
    '18obg4': 'arr-awq (Awanui)',
    '19obg5': 'arr-spq (Spring)',
    'greater_20obg6': 'arr-confq (Confined)'
}

# Convert to numeric
metadata_cols = ['iteration', 'total_runs', 'obs_realization', 'par_realization']
obs_groups = [col for col in phi_group.columns if col not in metadata_cols]
for col in obs_groups:
    phi_group[col] = pd.to_numeric(phi_group[col], errors='coerce')

# Calculate mean phi by iteration
phi_by_iter = phi_group.groupby('iteration')[obs_groups].mean()

# Read phi factors
phi_factors = pd.read_csv(os.path.join(master_dir, 'phi_factors.csv'), header=None, names=['group', 'weight'])

print("="*80)
print("ARR-H OBSERVATION GROUP ANALYSIS")
print("="*80)

# Detailed arr-h group performance
arr_h_data = []
for group, name in arr_h_groups.items():
    if group in phi_by_iter.columns:
        initial = phi_by_iter.loc[0, group]
        final = phi_by_iter.loc[phi_by_iter.index[-1], group]
        reduction = ((initial - final) / initial * 100) if initial > 0 else np.nan
        weight = phi_factors[phi_factors['group'] == group]['weight'].values[0] if group in phi_factors['group'].values else 0

        arr_h_data.append({
            'Group': group,
            'Name': name,
            'Weight': weight,
            'Initial_Phi': initial,
            'Final_Phi': final,
            'Reduction_%': reduction,
            'Ratio': final/initial if initial > 0 else np.nan
        })

arr_h_df = pd.DataFrame(arr_h_data)
print("\nARR-H GROUP SUMMARY:")
print(arr_h_df.to_string(index=False))

# Compare with non-arr-h groups
non_arr_h_groups = {
    '1obg0': 'Head obs 1',
    '2obg1': 'Head obs 2',
    '13obg2': 'Head obs 13',
    '3obg10': 'Recession rate 3',
    '22obg8': 'Recession rate 22',
    '23obg9': 'Recession rate 23',
    'greater_5obg12': 'Head obs greater_5',
    'less_6obg13': 'Head obs less_6',
    'less_16obg22': 'Head obs less_16'
}

print("\n" + "="*80)
print("COMPARISON WITH NON-ARR-H GROUPS:")
print("="*80)

comparison_data = []
for group, name in non_arr_h_groups.items():
    if group in phi_by_iter.columns:
        initial = phi_by_iter.loc[0, group]
        final = phi_by_iter.loc[phi_by_iter.index[-1], group]
        reduction = ((initial - final) / initial * 100) if initial > 0 else np.nan
        weight = phi_factors[phi_factors['group'] == group]['weight'].values[0] if group in phi_factors['group'].values else 0

        comparison_data.append({
            'Group': group,
            'Name': name,
            'Weight': weight,
            'Initial_Phi': initial,
            'Final_Phi': final,
            'Reduction_%': reduction,
            'Type': 'Non-ARR-H'
        })

for row in arr_h_data:
    comparison_data.append({
        'Group': row['Group'],
        'Name': row['Name'],
        'Weight': row['Weight'],
        'Initial_Phi': row['Initial_Phi'],
        'Final_Phi': row['Final_Phi'],
        'Reduction_%': row['Reduction_%'],
        'Type': 'ARR-H'
    })

comparison_df = pd.DataFrame(comparison_data).sort_values('Reduction_%', ascending=False)
print(comparison_df.to_string(index=False))

# Create detailed plots
fig = plt.figure(figsize=(16, 12))

# Plot 1: Phi evolution for arr-h groups
ax1 = plt.subplot(3, 2, 1)
for group, name in arr_h_groups.items():
    if group in phi_by_iter.columns:
        ax1.plot(phi_by_iter.index, phi_by_iter[group], marker='o', linewidth=2, label=name, markersize=6)
ax1.set_xlabel('Iteration', fontsize=11)
ax1.set_ylabel('Mean Phi', fontsize=11)
ax1.set_title('ARR-H Phi Evolution', fontsize=12, fontweight='bold')
ax1.legend(fontsize=9)
ax1.grid(True, alpha=0.3)
ax1.set_yscale('log')

# Plot 2: Phi evolution for selected non-arr-h groups
ax2 = plt.subplot(3, 2, 2)
selected_non_arr = ['1obg0', '2obg1', '13obg2', 'less_16obg22']
for group in selected_non_arr:
    if group in phi_by_iter.columns:
        name = non_arr_h_groups.get(group, group)
        ax2.plot(phi_by_iter.index, phi_by_iter[group], marker='s', linewidth=2, label=name, markersize=6)
ax2.set_xlabel('Iteration', fontsize=11)
ax2.set_ylabel('Mean Phi', fontsize=11)
ax2.set_title('Selected Head Obs Phi Evolution', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.set_yscale('log')

# Plot 3: Reduction percentage comparison
ax3 = plt.subplot(3, 2, 3)
colors = ['#e74c3c' if t == 'ARR-H' else '#3498db' for t in comparison_df['Type']]
bars = ax3.barh(range(len(comparison_df)), comparison_df['Reduction_%'], color=colors)
ax3.set_yticks(range(len(comparison_df)))
ax3.set_yticklabels(comparison_df['Name'], fontsize=9)
ax3.set_xlabel('Phi Reduction (%)', fontsize=11)
ax3.set_title('Phi Reduction from Initial', fontsize=12, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='x')
ax3.legend([plt.Rectangle((0,0),1,1, fc='#e74c3c'), plt.Rectangle((0,0),1,1, fc='#3498db')],
           ['ARR-H', 'Non-ARR-H'], loc='lower right', fontsize=9)

# Plot 4: Weight comparison
ax4 = plt.subplot(3, 2, 4)
arr_h_comp = comparison_df[comparison_df['Type'] == 'ARR-H']
non_arr_comp = comparison_df[comparison_df['Type'] == 'Non-ARR-H']
x = np.arange(3)
width = 0.35
ax4.bar(x - width/2, arr_h_comp.sort_values('Weight', ascending=False)['Weight'].values[:3],
        width, label='ARR-H', color='#e74c3c', alpha=0.8)
ax4.bar(x + width/2, non_arr_comp.sort_values('Weight', ascending=False)['Weight'].values[:3],
        width, label='Non-ARR-H', color='#3498db', alpha=0.8)
ax4.set_ylabel('Phi Weight', fontsize=11)
ax4.set_title('Observation Group Weights', fontsize=12, fontweight='bold')
ax4.set_xticks(x)
ax4.set_xticklabels(['Group 1', 'Group 2', 'Group 3'], fontsize=9)
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3, axis='y')

# Plot 5: Normalized phi (phi/initial_phi) for arr-h
ax5 = plt.subplot(3, 2, 5)
for group, name in arr_h_groups.items():
    if group in phi_by_iter.columns:
        initial = phi_by_iter.loc[0, group]
        normalized = phi_by_iter[group] / initial if initial > 0 else phi_by_iter[group]
        ax5.plot(phi_by_iter.index, normalized, marker='o', linewidth=2, label=name, markersize=6)
ax5.set_xlabel('Iteration', fontsize=11)
ax5.set_ylabel('Normalized Phi (Phi/Initial)', fontsize=11)
ax5.set_title('ARR-H Normalized Phi Evolution', fontsize=12, fontweight='bold')
ax5.legend(fontsize=9)
ax5.grid(True, alpha=0.3)
ax5.set_yscale('log')

# Plot 6: Final phi comparison
ax6 = plt.subplot(3, 2, 6)
final_comparison = comparison_df.sort_values('Final_Phi', ascending=True)
colors = ['#e74c3c' if t == 'ARR-H' else '#3498db' for t in final_comparison['Type']]
bars = ax6.barh(range(len(final_comparison)), final_comparison['Final_Phi'], color=colors)
ax6.set_yticks(range(len(final_comparison)))
ax6.set_yticklabels(final_comparison['Name'], fontsize=9)
ax6.set_xlabel('Final Phi (log scale)', fontsize=11)
ax6.set_title('Final Phi Values', fontsize=12, fontweight='bold')
ax6.set_xscale('log')
ax6.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('models/local_run33/arr_h_detailed_analysis.png', dpi=150, bbox_inches='tight')
print("\n" + "="*80)
print("Detailed figure saved to: models/local_run33/arr_h_detailed_analysis.png")
print("="*80)

# Calculate iteration-by-iteration changes
print("\n" + "="*80)
print("ITERATION-BY-ITERATION PHI CHANGES FOR ARR-H GROUPS:")
print("="*80)
for group, name in arr_h_groups.items():
    if group in phi_by_iter.columns:
        print(f"\n{name} ({group}):")
        print(f"{'Iteration':<12} {'Phi':<15} {'Change from prev':<20} {'% Change':<15}")
        print("-" * 65)
        for i, idx in enumerate(phi_by_iter.index):
            phi_val = phi_by_iter.loc[idx, group]
            if i == 0:
                print(f"{idx:<12} {phi_val:<15.2e} {'---':<20} {'---':<15}")
            else:
                prev_phi = phi_by_iter.loc[phi_by_iter.index[i-1], group]
                change = phi_val - prev_phi
                pct_change = (change / prev_phi * 100) if prev_phi > 0 else np.nan
                print(f"{idx:<12} {phi_val:<15.2e} {change:<20.2e} {pct_change:<15.2f}")
