"""
Analyze stuck realizations - why they don't improve during PEST++ IES
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

model_name = 'local_run33'
master_dir = f'models/{model_name}/pest/master_ies'

print("="*100)
print("ANALYZING STUCK REALIZATIONS")
print("="*100)

# Step 1: Identify stuck realizations from phi data
phi_group = pd.read_csv(f'{master_dir}/{model_name}.phi.group.csv')
metadata_cols = ['iteration', 'total_runs', 'obs_realization', 'par_realization']
obs_cols = [col for col in phi_group.columns if col not in metadata_cols]

# Calculate total phi
phi_totals = []
for idx, row in phi_group.iterrows():
    total = sum(pd.to_numeric(row[col], errors='coerce') for col in obs_cols
                if pd.notna(pd.to_numeric(row[col], errors='coerce')))
    phi_totals.append(total)
phi_group['total_phi'] = phi_totals

# Get initial and final phi
iter0 = phi_group[phi_group['iteration'] == 0][['obs_realization', 'total_phi']].copy()
iter0.columns = ['obs_realization', 'phi_iter0']

final_iter = phi_group['iteration'].max()
iter_final = phi_group[phi_group['iteration'] == final_iter][['obs_realization', 'total_phi']].copy()
iter_final.columns = ['obs_realization', 'phi_final']

comparison = iter0.merge(iter_final, on='obs_realization')
comparison['reduction_%'] = (comparison['phi_iter0'] - comparison['phi_final']) / comparison['phi_iter0'] * 100

# Identify stuck realizations (<50% reduction or negative reduction)
stuck = comparison[comparison['reduction_%'] < 50].copy()
stuck = stuck.sort_values('reduction_%')

print(f"\nFound {len(stuck)} stuck realizations (out of {len(comparison)})")
print(f"That's {len(stuck)/len(comparison)*100:.1f}% of all realizations!")

# Get worst 10
worst_10 = stuck.head(10)
print(f"\n10 Worst Performing Realizations:")
print(worst_10[['obs_realization', 'phi_iter0', 'phi_final', 'reduction_%']].to_string(index=False))

# Step 2: Read parameter data
print(f"\n{'='*100}")
print("READING PARAMETER DATA...")
print("="*100)

par0 = pd.read_csv(f'{master_dir}/{model_name}.0.par.csv', low_memory=False)
par_final = pd.read_csv(f'{master_dir}/{model_name}.{final_iter}.par.csv', low_memory=False)

print(f"Iteration 0: {par0.shape[0]} realizations, {par0.shape[1]-1} parameters")
print(f"Iteration {final_iter}: {par_final.shape[0]} realizations, {par_final.shape[1]-1} parameters")

# Get parameter columns
param_cols = [col for col in par0.columns if col != 'real_name']

# Step 3: Calculate parameter changes for stuck vs good realizations
print(f"\n{'='*100}")
print("COMPARING PARAMETER CHANGES: STUCK vs GOOD REALIZATIONS")
print("="*100)

# Define good realizations (>95% reduction)
good = comparison[comparison['reduction_%'] > 95].copy()
print(f"Good realizations (>95% reduction): {len(good)}")

# Get parameter values for stuck and good realizations
stuck_reals = stuck['obs_realization'].values
good_reals = good['obs_realization'].head(20).values  # Sample 20 good ones

# For stuck realizations, calculate parameter change magnitude
def calc_param_change(real_id, par0, par_final, param_cols):
    """Calculate parameter change for a single realization"""
    try:
        p0 = par0[par0['real_name'] == real_id][param_cols].values[0]
        pf = par_final[par_final['real_name'] == real_id][param_cols].values[0]

        # Calculate various metrics
        abs_change = np.abs(pf - p0)
        rel_change = np.abs((pf - p0) / (p0 + 1e-10))  # Relative change

        return {
            'mean_abs_change': np.mean(abs_change),
            'max_abs_change': np.max(abs_change),
            'mean_rel_change': np.mean(rel_change),
            'max_rel_change': np.max(rel_change),
            'num_changed': np.sum(abs_change > 0.01),  # Parameters that changed >0.01
            'pct_changed': np.sum(abs_change > 0.01) / len(param_cols) * 100
        }
    except:
        return None

print("\nAnalyzing parameter changes for worst 10 stuck realizations...")
stuck_param_changes = []
for real_id in worst_10['obs_realization'].values:
    change = calc_param_change(real_id, par0, par_final, param_cols)
    if change:
        change['real_id'] = real_id
        change['type'] = 'stuck'
        stuck_param_changes.append(change)

print("Analyzing parameter changes for 20 good realizations...")
good_param_changes = []
for real_id in good_reals:
    change = calc_param_change(real_id, par0, par_final, param_cols)
    if change:
        change['real_id'] = real_id
        change['type'] = 'good'
        good_param_changes.append(change)

stuck_df = pd.DataFrame(stuck_param_changes)
good_df = pd.DataFrame(good_param_changes)

print(f"\n{'='*100}")
print("PARAMETER CHANGE STATISTICS")
print("="*100)

print("\nSTUCK REALIZATIONS (worst 10):")
print(stuck_df[['real_id', 'mean_abs_change', 'max_abs_change', 'pct_changed']].to_string(index=False))

print("\nGOOD REALIZATIONS (sample of 20):")
print(good_df[['real_id', 'mean_abs_change', 'max_abs_change', 'pct_changed']].to_string(index=False))

print(f"\n{'='*100}")
print("COMPARISON SUMMARY")
print("="*100)

summary = pd.DataFrame({
    'Metric': ['Mean Abs Change', 'Max Abs Change', 'Mean Rel Change', 'Max Rel Change', '% Params Changed'],
    'Stuck (worst 10)': [
        stuck_df['mean_abs_change'].mean(),
        stuck_df['max_abs_change'].mean(),
        stuck_df['mean_rel_change'].mean(),
        stuck_df['max_rel_change'].mean(),
        stuck_df['pct_changed'].mean()
    ],
    'Good (sample 20)': [
        good_df['mean_abs_change'].mean(),
        good_df['max_abs_change'].mean(),
        good_df['mean_rel_change'].mean(),
        good_df['max_rel_change'].mean(),
        good_df['pct_changed'].mean()
    ]
})
summary['Ratio (Stuck/Good)'] = summary['Stuck (worst 10)'] / summary['Good (sample 20)']
print(summary.to_string(index=False))

# Step 4: Check if stuck realizations hit parameter bounds
print(f"\n{'='*100}")
print("CHECKING PARAMETER BOUNDS")
print("="*100)

# Read parameter bounds
par_data = pd.read_csv(f'{master_dir}/{model_name}.par_data.csv')
param_bounds = dict(zip(par_data['parnme'], zip(par_data['parlbnd'], par_data['parubnd'])))

def check_bounds(real_id, par_final, param_cols, param_bounds):
    """Check how many parameters are at bounds"""
    try:
        pf = par_final[par_final['real_name'] == real_id][param_cols].values[0]

        at_lower = 0
        at_upper = 0
        near_lower = 0
        near_upper = 0

        for i, col in enumerate(param_cols):
            lb, ub = param_bounds[col]
            val = pf[i]

            if np.abs(val - lb) < 1e-6:
                at_lower += 1
            elif np.abs(val - ub) < 1e-6:
                at_upper += 1
            elif np.abs(val - lb) < 0.1 * (ub - lb):
                near_lower += 1
            elif np.abs(val - ub) < 0.1 * (ub - lb):
                near_upper += 1

        return {
            'at_lower': at_lower,
            'at_upper': at_upper,
            'near_lower': near_lower,
            'near_upper': near_upper,
            'pct_at_bounds': (at_lower + at_upper) / len(param_cols) * 100
        }
    except:
        return None

print("\nChecking bounds for worst 10 stuck realizations...")
stuck_bounds = []
for real_id in worst_10['obs_realization'].values:
    bounds = check_bounds(real_id, par_final, param_cols, param_bounds)
    if bounds:
        bounds['real_id'] = real_id
        stuck_bounds.append(bounds)

stuck_bounds_df = pd.DataFrame(stuck_bounds)
print(stuck_bounds_df[['real_id', 'at_lower', 'at_upper', 'pct_at_bounds']].to_string(index=False))

print(f"\n{'='*100}")
print("CONCLUSIONS")
print("="*100)

print(f"""
1. Stuck Realizations: {len(stuck)} out of {len(comparison)} ({len(stuck)/len(comparison)*100:.1f}%)

2. Parameter Changes:
   - Stuck realizations change parameters by: {stuck_df['mean_abs_change'].mean():.4f} (mean)
   - Good realizations change parameters by: {good_df['mean_abs_change'].mean():.4f} (mean)
   - Ratio: {stuck_df['mean_abs_change'].mean() / good_df['mean_abs_change'].mean():.2f}x

3. Percentage of Parameters Changed:
   - Stuck: {stuck_df['pct_changed'].mean():.1f}%
   - Good: {good_df['pct_changed'].mean():.1f}%

4. Parameters at Bounds:
   - Mean % at bounds for stuck: {stuck_bounds_df['pct_at_bounds'].mean():.2f}%

DIAGNOSIS:
""")

if stuck_df['mean_abs_change'].mean() < 0.1 * good_df['mean_abs_change'].mean():
    print("⚠️  STUCK REALIZATIONS ARE NOT BEING UPDATED!")
    print("   Parameters are barely changing between iterations.")
    print("   This suggests:")
    print("   - Low ensemble variance (realizations too similar)")
    print("   - Insufficient localization")
    print("   - Parameters may have converged to local minimum")

if stuck_bounds_df['pct_at_bounds'].mean() > 10:
    print("⚠️  MANY PARAMETERS AT BOUNDS!")
    print(f"   {stuck_bounds_df['pct_at_bounds'].mean():.1f}% of parameters are at bounds")
    print("   This suggests:")
    print("   - Parameter bounds may be too restrictive")
    print("   - True parameter values may be outside allowed range")

print("\nSaving results...")
stuck.to_csv(f'models/{model_name}/stuck_realizations.csv', index=False)
comparison.to_csv(f'models/{model_name}/realization_performance.csv', index=False)

print(f"\nFiles saved:")
print(f"  - models/{model_name}/stuck_realizations.csv")
print(f"  - models/{model_name}/realization_performance.csv")
