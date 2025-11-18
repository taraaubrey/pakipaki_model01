"""
Final comprehensive summary of arr-h observations
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read data
master_dir = 'models/local_run33/pest/master_ies'
phi_group = pd.read_csv(f'{master_dir}/local_run33.phi.group.csv')
phi_actual = pd.read_csv(f'{master_dir}/local_run33.phi.actual.csv')

# Convert to numeric
metadata_cols = ['iteration', 'total_runs', 'obs_realization', 'par_realization']
obs_groups = [col for col in phi_group.columns if col not in metadata_cols]
for col in obs_groups:
    phi_group[col] = pd.to_numeric(phi_group[col], errors='coerce')

phi_by_iter = phi_group.groupby('iteration')[obs_groups].mean()

# Define groups
arr_h_groups = {
    '18obg4': 'Awanui Stream',
    '19obg5': 'Springs',
    'greater_20obg6': 'Confined Area'
}

head_groups = {
    '1obg0': 'Head-1',
    '2obg1': 'Head-2',
    '13obg2': 'Head-13',
    'less_16obg22': 'Head-less16'
}

recession_groups = {
    '3obg10': 'Recession-3',
    '22obg8': 'Recession-22',
    '23obg9': 'Recession-23'
}

# Create comprehensive figure
fig = plt.figure(figsize=(18, 14))
gs = fig.add_gridspec(4, 3, hspace=0.3, wspace=0.3)

# Plot 1: Total phi evolution
ax1 = fig.add_subplot(gs[0, :])
ax1.plot(phi_actual['iteration'], phi_actual['mean'], 'o-', linewidth=3, markersize=8, color='black', label='Mean Total Phi')
ax1.fill_between(phi_actual['iteration'],
                  phi_actual['mean'] - phi_actual['standard_deviation'],
                  phi_actual['mean'] + phi_actual['standard_deviation'],
                  alpha=0.2, color='black')
ax1.axvline(x=4, color='red', linestyle='--', linewidth=2, alpha=0.5, label='Iteration 4 (spike)')
ax1.set_xlabel('Iteration', fontsize=12, fontweight='bold')
ax1.set_ylabel('Total Phi', fontsize=12, fontweight='bold')
ax1.set_title('Total Phi Evolution Across All Observations', fontsize=14, fontweight='bold')
ax1.set_yscale('log')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Plot 2: ARR-H phi evolution
ax2 = fig.add_subplot(gs[1, 0])
for group, name in arr_h_groups.items():
    if group in phi_by_iter.columns:
        ax2.plot(phi_by_iter.index, phi_by_iter[group], marker='o', linewidth=2, label=name, markersize=5)
ax2.axvline(x=4, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
ax2.set_xlabel('Iteration', fontsize=11)
ax2.set_ylabel('Mean Phi', fontsize=11)
ax2.set_title('ARR-H Observations', fontsize=12, fontweight='bold')
ax2.legend(fontsize=9)
ax2.grid(True, alpha=0.3)
ax2.set_yscale('log')

# Plot 3: Head obs phi evolution
ax3 = fig.add_subplot(gs[1, 1])
for group, name in head_groups.items():
    if group in phi_by_iter.columns:
        ax3.plot(phi_by_iter.index, phi_by_iter[group], marker='s', linewidth=2, label=name, markersize=5)
ax3.axvline(x=4, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
ax3.set_xlabel('Iteration', fontsize=11)
ax3.set_ylabel('Mean Phi', fontsize=11)
ax3.set_title('Head Observations', fontsize=12, fontweight='bold')
ax3.legend(fontsize=9)
ax3.grid(True, alpha=0.3)
ax3.set_yscale('log')

# Plot 4: Recession obs phi evolution
ax4 = fig.add_subplot(gs[1, 2])
for group, name in recession_groups.items():
    if group in phi_by_iter.columns:
        ax4.plot(phi_by_iter.index, phi_by_iter[group], marker='^', linewidth=2, label=name, markersize=5)
ax4.axvline(x=4, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
ax4.set_xlabel('Iteration', fontsize=11)
ax4.set_ylabel('Mean Phi', fontsize=11)
ax4.set_title('Recession Rate Observations', fontsize=12, fontweight='bold')
ax4.legend(fontsize=9)
ax4.grid(True, alpha=0.3)
ax4.set_yscale('log')

# Plot 5: Normalized phi for all groups
ax5 = fig.add_subplot(gs[2, :2])
all_groups = {**arr_h_groups, **head_groups, **recession_groups}
colors_map = {**{k: '#e74c3c' for k in arr_h_groups},
              **{k: '#3498db' for k in head_groups},
              **{k: '#2ecc71' for k in recession_groups}}
for group, name in all_groups.items():
    if group in phi_by_iter.columns:
        initial = phi_by_iter.loc[0, group]
        normalized = phi_by_iter[group] / initial if initial > 0 else phi_by_iter[group]
        ax5.plot(phi_by_iter.index, normalized, linewidth=2, label=name, color=colors_map[group], alpha=0.7)
ax5.axvline(x=4, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
ax5.axhline(y=1, color='gray', linestyle=':', linewidth=1)
ax5.set_xlabel('Iteration', fontsize=11)
ax5.set_ylabel('Normalized Phi (Phi/Initial)', fontsize=11)
ax5.set_title('Normalized Phi Evolution (All Groups)', fontsize=12, fontweight='bold')
ax5.legend(ncol=2, fontsize=8, loc='upper right')
ax5.grid(True, alpha=0.3)
ax5.set_yscale('log')

# Plot 6: Final phi comparison
ax6 = fig.add_subplot(gs[2, 2])
final_phi_data = []
for group, name in all_groups.items():
    if group in phi_by_iter.columns:
        group_type = 'ARR-H' if group in arr_h_groups else ('Head' if group in head_groups else 'Recession')
        final_phi_data.append({
            'name': name,
            'phi': phi_by_iter.loc[phi_by_iter.index[-1], group],
            'type': group_type
        })
final_df = pd.DataFrame(final_phi_data).sort_values('phi', ascending=True)
colors = [colors_map[k] for k, v in all_groups.items() if v in final_df['name'].values]
bars = ax6.barh(range(len(final_df)), final_df['phi'],
                color=['#e74c3c' if t == 'ARR-H' else '#3498db' if t == 'Head' else '#2ecc71' for t in final_df['type']])
ax6.set_yticks(range(len(final_df)))
ax6.set_yticklabels(final_df['name'], fontsize=9)
ax6.set_xlabel('Final Phi (log scale)', fontsize=11)
ax6.set_title('Final Phi Values', fontsize=12, fontweight='bold')
ax6.set_xscale('log')
ax6.grid(True, alpha=0.3, axis='x')

# Plot 7: Phi reduction percentages
ax7 = fig.add_subplot(gs[3, :])
reduction_data = []
for group, name in all_groups.items():
    if group in phi_by_iter.columns:
        initial = phi_by_iter.loc[0, group]
        final = phi_by_iter.loc[phi_by_iter.index[-1], group]
        reduction = ((initial - final) / initial * 100) if initial > 0 else 0
        group_type = 'ARR-H' if group in arr_h_groups else ('Head' if group in head_groups else 'Recession')
        reduction_data.append({
            'name': name,
            'reduction': reduction,
            'type': group_type
        })
reduction_df = pd.DataFrame(reduction_data).sort_values('reduction', ascending=False)
colors_bar = ['#e74c3c' if t == 'ARR-H' else '#3498db' if t == 'Head' else '#2ecc71' for t in reduction_df['type']]
bars = ax7.bar(range(len(reduction_df)), reduction_df['reduction'], color=colors_bar, alpha=0.8, edgecolor='black')
ax7.set_xticks(range(len(reduction_df)))
ax7.set_xticklabels(reduction_df['name'], rotation=45, ha='right', fontsize=9)
ax7.set_ylabel('Phi Reduction (%)', fontsize=11)
ax7.set_title('Phi Reduction from Initial by Observation Group', fontsize=12, fontweight='bold')
ax7.grid(True, alpha=0.3, axis='y')
ax7.axhline(y=99, color='gray', linestyle='--', linewidth=1, alpha=0.5)
ax7.axhline(y=90, color='gray', linestyle=':', linewidth=1, alpha=0.5)
# Add legend
from matplotlib.patches import Rectangle
legend_elements = [Rectangle((0,0),1,1, fc='#e74c3c', label='ARR-H'),
                   Rectangle((0,0),1,1, fc='#3498db', label='Head Obs'),
                   Rectangle((0,0),1,1, fc='#2ecc71', label='Recession')]
ax7.legend(handles=legend_elements, loc='lower left', fontsize=10)

plt.savefig('models/local_run33/COMPREHENSIVE_ANALYSIS.png', dpi=200, bbox_inches='tight')
print('Comprehensive figure saved to: models/local_run33/COMPREHENSIVE_ANALYSIS.png')

# Print summary table
print('\n' + '='*100)
print('FINAL SUMMARY TABLE')
print('='*100)
summary_data = []
for group, name in all_groups.items():
    if group in phi_by_iter.columns:
        initial = phi_by_iter.loc[0, group]
        final = phi_by_iter.loc[phi_by_iter.index[-1], group]
        reduction = ((initial - final) / initial * 100) if initial > 0 else 0
        group_type = 'ARR-H' if group in arr_h_groups else ('Head' if group in head_groups else 'Recession')
        summary_data.append({
            'Type': group_type,
            'Name': name,
            'Initial_Phi': f'{initial:.2e}',
            'Final_Phi': f'{final:.2e}',
            'Reduction_%': f'{reduction:.2f}'
        })
summary_table = pd.DataFrame(summary_data)
print(summary_table.to_string(index=False))
