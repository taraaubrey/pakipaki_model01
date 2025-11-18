"""
Analyze PEST++ IES phi results for arr-h and other observation groups.
Run this after e_pest_IES_HM.py to understand observation group performance.

Usage:
    python scripts/f_analyze_phi.py [model_name]

    If model_name is not provided, uses MODEL_NAME from setup.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

# Import setup to get model name
from setup import MODEL_NAME as DEFAULT_MODEL_NAME

def get_base_obgnme(obgnme):
    """Remove leading numbers from observation group names like '16obg4' -> 'obg4'"""
    match = re.match(r'^(\d+|greater_\d+|less_\d+)(.*)$', str(obgnme))
    if match and match.group(2):
        return match.group(2)
    return obgnme


def analyze_phi(model_name):
    """
    Analyze phi results for a given model run.

    Args:
        model_name: Name of the model run (e.g., 'local_run33')
    """
    print("="*100)
    print(f"PEST++ IES PHI ANALYSIS: {model_name}")
    print("="*100)

    # Set paths
    master_dir = os.path.join('models', model_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', model_name)

    if not os.path.exists(master_dir):
        print(f"ERROR: Directory not found: {master_dir}")
        print("Make sure you've run e_pest_IES_HM.py first!")
        return

    # Read phi data
    phi_group_file = os.path.join(master_dir, f'{model_name}.phi.group.csv')
    phi_actual_file = os.path.join(master_dir, f'{model_name}.phi.actual.csv')
    phi_factors_file = os.path.join(master_dir, 'phi_factors.csv')
    obs_data_file = os.path.join(master_dir, f'{model_name}.obs_data.csv')

    print(f"\nReading phi data from: {master_dir}")
    phi_group = pd.read_csv(phi_group_file)
    phi_actual = pd.read_csv(phi_actual_file)
    phi_factors = pd.read_csv(phi_factors_file, header=None, names=['group', 'weight'])

    # Read obs_data to map observation types
    # Read just the columns we need to save memory
    print(f"Reading observation data...")
    obs_data = pd.read_csv(obs_data_file, usecols=['obgnme', 'oname', 'usecol'], low_memory=False)

    # Create mapping of obgnme to observation type (oname + usecol)
    obs_mapping = obs_data[['obgnme', 'oname', 'usecol']].drop_duplicates()
    # Combine oname and usecol for better description
    obs_mapping['full_name'] = obs_mapping.apply(
        lambda x: f"{x['oname']}:{x['usecol']}" if pd.notna(x['usecol']) and str(x['usecol']) != '0'
        else x['oname'], axis=1
    )
    obs_type_map = dict(zip(obs_mapping['obgnme'], obs_mapping['full_name']))

    # Identify observation groups by type
    arr_h_groups = {}
    head_groups = {}
    recession_groups = {}
    flux_groups = {}
    budget_groups = {}
    other_groups = {}

    for obgnme, oname in obs_type_map.items():
        oname_lower = str(oname).lower()

        if 'arr' in oname_lower:
            arr_h_groups[obgnme] = oname
        elif 'head' in oname_lower or 'ts-heads' in oname_lower:
            head_groups[obgnme] = oname
        elif 'rec' in oname_lower or 'recession' in oname_lower:
            # Exclude budget:recharge from recession group
            if 'budget' not in oname_lower:
                recession_groups[obgnme] = oname
            else:
                budget_groups[obgnme] = oname
        elif 'flux' in oname_lower or 'flow' in oname_lower or 'sprflux' in oname_lower:
            flux_groups[obgnme] = oname
        elif 'budget' in oname_lower:
            budget_groups[obgnme] = oname
        else:
            other_groups[obgnme] = oname

    print(f"\nIdentified observation groups:")
    print(f"\n  ARR-H groups ({len(arr_h_groups)}):")
    for group, oname in sorted(arr_h_groups.items()):
        weight = phi_factors[phi_factors['group'] == group]['weight'].values
        weight_val = weight[0] if len(weight) > 0 else 0.0
        print(f"    {group}: {oname} (weight={weight_val})")

    print(f"\n  Head groups ({len(head_groups)}):")
    for group, oname in sorted(head_groups.items()):
        weight = phi_factors[phi_factors['group'] == group]['weight'].values
        weight_val = weight[0] if len(weight) > 0 else 0.0
        print(f"    {group}: {oname} (weight={weight_val})")

    print(f"\n  Recession groups ({len(recession_groups)}):")
    for group, oname in sorted(recession_groups.items()):
        weight = phi_factors[phi_factors['group'] == group]['weight'].values
        weight_val = weight[0] if len(weight) > 0 else 0.0
        print(f"    {group}: {oname} (weight={weight_val})")

    if len(flux_groups) > 0:
        print(f"\n  Flux groups ({len(flux_groups)}):")
        for group, oname in sorted(flux_groups.items()):
            weight = phi_factors[phi_factors['group'] == group]['weight'].values
            weight_val = weight[0] if len(weight) > 0 else 0.0
            print(f"    {group}: {oname} (weight={weight_val})")

    if len(budget_groups) > 0:
        print(f"\n  Budget groups ({len(budget_groups)}):")
        for group, oname in sorted(budget_groups.items()):
            weight = phi_factors[phi_factors['group'] == group]['weight'].values
            weight_val = weight[0] if len(weight) > 0 else 0.0
            print(f"    {group}: {oname} (weight={weight_val})")

    # Process phi_group data
    metadata_cols = ['iteration', 'total_runs', 'obs_realization', 'par_realization']
    obs_group_cols = [col for col in phi_group.columns if col not in metadata_cols]

    # Convert to numeric
    for col in obs_group_cols:
        phi_group[col] = pd.to_numeric(phi_group[col], errors='coerce')

    # Calculate mean phi by iteration
    phi_by_iter = phi_group.groupby('iteration')[obs_group_cols].mean()

    # Calculate statistics for each group type
    def get_group_stats(groups_dict, group_type_name):
        """Calculate stats for a group type"""
        stats = []
        for group in groups_dict.keys():
            if group not in phi_by_iter.columns:
                continue

            initial = phi_by_iter.loc[0, group]
            final = phi_by_iter.loc[phi_by_iter.index[-1], group]
            reduction = ((initial - final) / initial * 100) if initial > 0 else np.nan
            weight = phi_factors[phi_factors['group'] == group]['weight'].values
            weight_val = weight[0] if len(weight) > 0 else 0.0

            stats.append({
                'Type': group_type_name,
                'Group': group,
                'Obs_Name': groups_dict[group],
                'Weight': weight_val,
                'Initial_Phi': initial,
                'Final_Phi': final,
                'Reduction_%': reduction
            })
        return stats

    # Compile all statistics
    all_stats = []
    all_stats.extend(get_group_stats(arr_h_groups, 'ARR-H'))
    all_stats.extend(get_group_stats(head_groups, 'Head'))
    all_stats.extend(get_group_stats(recession_groups, 'Recession'))
    all_stats.extend(get_group_stats(flux_groups, 'Flux'))
    all_stats.extend(get_group_stats(budget_groups, 'Budget'))
    all_stats.extend(get_group_stats(other_groups, 'Other'))

    stats_df = pd.DataFrame(all_stats).sort_values('Reduction_%', ascending=False)

    print("\n" + "="*100)
    print("OBSERVATION GROUP PERFORMANCE SUMMARY")
    print("="*100)
    print(stats_df.to_string(index=False))

    # Calculate final contribution to total phi
    final_iter = phi_by_iter.index[-1]
    final_total = phi_by_iter.loc[final_iter].sum()
    stats_df['Final_Contribution_%'] = (stats_df['Final_Phi'] / final_total * 100)

    print("\n" + "="*100)
    print("FINAL PHI CONTRIBUTION (% of total residual)")
    print("="*100)
    for obs_type in ['ARR-H', 'Head', 'Recession', 'Flux', 'Budget', 'Other']:
        contrib = stats_df[stats_df['Type'] == obs_type]['Final_Contribution_%'].sum()
        if contrib > 0:
            print(f"{obs_type:20s}: {contrib:6.2f}%")

    # Total phi evolution
    print("\n" + "="*100)
    print("TOTAL PHI EVOLUTION")
    print("="*100)
    phi_actual['reduction_%'] = (phi_actual.loc[0, 'mean'] - phi_actual['mean']) / phi_actual.loc[0, 'mean'] * 100
    print(phi_actual[['iteration', 'mean', 'standard_deviation', 'reduction_%']].to_string(index=False))

    # Create comprehensive visualization
    print("\n" + "="*100)
    print("CREATING VISUALIZATIONS...")
    print("="*100)

    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(4, 3, hspace=0.35, wspace=0.35)

    # Plot 1: Total phi evolution
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(phi_actual['iteration'], phi_actual['mean'], 'o-', linewidth=3, markersize=8,
             color='black', label='Mean Total Phi', zorder=3)
    ax1.fill_between(phi_actual['iteration'],
                      phi_actual['mean'] - phi_actual['standard_deviation'],
                      phi_actual['mean'] + phi_actual['standard_deviation'],
                      alpha=0.2, color='black', label='±1 Std Dev')

    # Mark reinflation iterations
    reinflate_iters = phi_actual[phi_actual['reduction_%'] < 0]['iteration'].values
    for iter_num in reinflate_iters:
        ax1.axvline(x=iter_num, color='red', linestyle='--', linewidth=1.5, alpha=0.5)
    if len(reinflate_iters) > 0:
        ax1.plot([], [], 'r--', label='Reinflation', linewidth=1.5, alpha=0.5)

    ax1.set_xlabel('Iteration', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Total Phi', fontsize=12, fontweight='bold')
    ax1.set_title(f'Total Phi Evolution: {model_name}', fontsize=14, fontweight='bold')
    ax1.set_yscale('log')
    ax1.legend(fontsize=10, loc='upper right')
    ax1.grid(True, alpha=0.3)

    # Define color scheme function (same as used in plots 6 & 7)
    def get_obs_color(obs_name):
        """Assign colors based on observation subcategory"""
        obs_lower = obs_name.lower()
        # ARR stream/spring (awq, spq, confq) - shades of red/pink
        if 'arr-awq' in obs_lower or 'arr-spq' in obs_lower or 'arr-confq' in obs_lower:
            if 'awq' in obs_lower:
                return '#e74c3c'  # Red for Awanui
            elif 'spq' in obs_lower:
                return '#c0392b'  # Dark red for springs
            else:
                return '#e67e22'  # Orange for confined arr
        # ARR-H and ts-heads - shades of blue
        elif 'arr-h' in obs_lower or 'ts-heads' in obs_lower:
            if 'arr-h' in obs_lower:
                return '#3498db'  # Blue for arr-h
            else:
                return '#2980b9'  # Dark blue for ts-heads
        # Recession - shades of green
        elif 'recession' in obs_lower:
            return '#2ecc71'  # Green for recession
        # Budget - shades of purple/gray
        elif 'budget' in obs_lower:
            return '#9b59b6'  # Purple for budget
        else:
            return '#95a5a6'  # Gray for other

    # Plots 2-4: Phi evolution by group type (only non-zero weights)
    group_configs = [
        (arr_h_groups, 'ARR-H Observations', 'o', gs[1, 0]),
        (head_groups, 'Head Observations', 's', gs[1, 1]),
        (recession_groups, 'Recession Observations', '^', gs[1, 2])
    ]

    for groups, title, marker, subplot_loc in group_configs:
        ax = fig.add_subplot(subplot_loc)
        if len(groups) > 0:
            for group, oname in groups.items():
                if group in phi_by_iter.columns:
                    # Only plot groups with non-zero weights
                    weight = phi_factors[phi_factors['group'] == group]['weight'].values
                    weight_val = weight[0] if len(weight) > 0 else 0.0
                    if weight_val > 0:
                        color = get_obs_color(oname)
                        ax.plot(phi_by_iter.index, phi_by_iter[group], marker=marker,
                               linewidth=2, label=f"{oname}", markersize=5, alpha=0.8, color=color)

            for iter_num in reinflate_iters:
                ax.axvline(x=iter_num, color='red', linestyle='--', linewidth=1, alpha=0.3)

            ax.legend(fontsize=8, loc='best')
        else:
            ax.text(0.5, 0.5, 'No observations', ha='center', va='center',
                   transform=ax.transAxes, fontsize=12, color='gray')

        ax.set_xlabel('Iteration', fontsize=11)
        ax.set_ylabel('Mean Phi', fontsize=11)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.set_yscale('log')

    # Plot 5: Normalized phi comparison (only non-zero weights)
    ax5 = fig.add_subplot(gs[2, :2])
    all_groups = {**arr_h_groups, **head_groups, **recession_groups}

    for group, oname in all_groups.items():
        if group in phi_by_iter.columns:
            # Only plot groups with non-zero weights
            weight = phi_factors[phi_factors['group'] == group]['weight'].values
            weight_val = weight[0] if len(weight) > 0 else 0.0
            if weight_val > 0:
                initial = phi_by_iter.loc[0, group]
                normalized = phi_by_iter[group] / initial if initial > 0 else phi_by_iter[group]
                color = get_obs_color(oname)
                ax5.plot(phi_by_iter.index, normalized, linewidth=2, label=oname,
                        color=color, alpha=0.7)

    for iter_num in reinflate_iters:
        ax5.axvline(x=iter_num, color='red', linestyle='--', linewidth=1, alpha=0.3)

    ax5.axhline(y=1, color='gray', linestyle=':', linewidth=1)
    ax5.set_xlabel('Iteration', fontsize=11)
    ax5.set_ylabel('Normalized Phi (Phi/Initial)', fontsize=11)
    ax5.set_title('Normalized Phi Evolution (All Groups)', fontsize=12, fontweight='bold')
    ax5.legend(ncol=2, fontsize=8, loc='upper right')
    ax5.grid(True, alpha=0.3)
    ax5.set_yscale('log')

    # Plot 6: Final phi comparison (only non-zero weighted obs)
    ax6 = fig.add_subplot(gs[2, 2])
    # Filter to only show observations with non-zero weights
    final_df = stats_df[stats_df['Weight'] > 0].sort_values('Final_Phi', ascending=True)

    if len(final_df) > 0:
        bar_colors = [get_obs_color(name) for name in final_df['Obs_Name']]
        ax6.barh(range(len(final_df)), final_df['Final_Phi'], color=bar_colors, alpha=0.8, edgecolor='black')
        ax6.set_yticks(range(len(final_df)))
        ax6.set_yticklabels(final_df['Obs_Name'], fontsize=8)
        ax6.set_xlabel('Final Phi (log scale)', fontsize=11)
        ax6.set_title('Final Phi Values (Non-Zero Weights Only)', fontsize=12, fontweight='bold')
        ax6.set_xscale('log')
        ax6.grid(True, alpha=0.3, axis='x')
    else:
        ax6.text(0.5, 0.5, 'No non-zero weighted observations', ha='center', va='center',
                transform=ax6.transAxes, fontsize=12, color='gray')

    # Plot 7: Phi reduction percentages (only non-zero weighted obs)
    ax7 = fig.add_subplot(gs[3, :])
    # Filter to only show observations with non-zero weights
    reduction_df = stats_df[stats_df['Weight'] > 0].sort_values('Reduction_%', ascending=False)

    if len(reduction_df) > 0:
        # Use the same color scheme as Plot 6
        bar_colors = [get_obs_color(name) for name in reduction_df['Obs_Name']]

        bars = ax7.bar(range(len(reduction_df)), reduction_df['Reduction_%'],
                       color=bar_colors, alpha=0.8, edgecolor='black')
        ax7.set_xticks(range(len(reduction_df)))
        ax7.set_xticklabels(reduction_df['Obs_Name'], rotation=45, ha='right', fontsize=9)
        ax7.set_ylabel('Phi Reduction (%)', fontsize=11)
        ax7.set_title('Phi Reduction from Initial by Observation Group (Non-Zero Weights Only)',
                     fontsize=12, fontweight='bold')
        ax7.grid(True, alpha=0.3, axis='y')
        ax7.axhline(y=99, color='gray', linestyle='--', linewidth=1, alpha=0.5, label='99%')
        ax7.axhline(y=90, color='gray', linestyle=':', linewidth=1, alpha=0.5, label='90%')

        # Add custom legend for color groups
        from matplotlib.patches import Rectangle
        legend_elements = [
            Rectangle((0,0),1,1, fc='#e74c3c', ec='black', label='arr-awq'),
            Rectangle((0,0),1,1, fc='#c0392b', ec='black', label='arr-spq'),
            Rectangle((0,0),1,1, fc='#e67e22', ec='black', label='arr-confq'),
            Rectangle((0,0),1,1, fc='#3498db', ec='black', label='arr-h'),
            Rectangle((0,0),1,1, fc='#2980b9', ec='black', label='ts-heads'),
            Rectangle((0,0),1,1, fc='#2ecc71', ec='black', label='recession'),
            Rectangle((0,0),1,1, fc='#9b59b6', ec='black', label='budget'),
        ]
        # Only show legend items that are actually present
        present_colors = set(bar_colors)
        legend_map = {
            '#e74c3c': 'arr-awq', '#c0392b': 'arr-spq', '#e67e22': 'arr-confq',
            '#3498db': 'arr-h', '#2980b9': 'ts-heads', '#2ecc71': 'recession',
            '#9b59b6': 'budget'
        }
        filtered_legend = [Rectangle((0,0),1,1, fc=c, ec='black', label=legend_map[c])
                          for c in present_colors if c in legend_map]
        if filtered_legend:
            ax7.legend(handles=filtered_legend, loc='lower left', fontsize=9, ncol=3)
    else:
        ax7.text(0.5, 0.5, 'No non-zero weighted observations', ha='center', va='center',
                transform=ax7.transAxes, fontsize=12, color='gray')

    # Save figure
    fig_path = os.path.join(output_dir, f'{model_name}_phi_analysis.png')
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"Figure saved to: {fig_path}")

    # Save summary to CSV
    csv_path = os.path.join(output_dir, f'{model_name}_phi_summary.csv')
    stats_df.to_csv(csv_path, index=False)
    print(f"Summary table saved to: {csv_path}")

    # Create markdown summary
    md_path = os.path.join(output_dir, f'{model_name}_PHI_ANALYSIS.md')
    with open(md_path, 'w') as f:
        f.write(f"# Phi Analysis Summary: {model_name}\n\n")
        f.write("## Overview\n\n")
        f.write(f"Total iterations completed: {phi_actual['iteration'].max()}\n\n")
        f.write(f"Initial phi (mean): {phi_actual.loc[0, 'mean']:.2e}\n\n")
        f.write(f"Final phi (mean): {phi_actual.loc[phi_actual.index[-1], 'mean']:.2e}\n\n")
        f.write(f"Overall reduction: {phi_actual.loc[phi_actual.index[-1], 'reduction_%']:.2f}%\n\n")

        if len(reinflate_iters) > 0:
            f.write(f"Reinflation iterations: {', '.join(map(str, reinflate_iters))}\n\n")

        f.write("## ARR-H Observation Performance\n\n")
        arr_h_stats = stats_df[stats_df['Type'] == 'ARR-H']
        if len(arr_h_stats) > 0:
            f.write(arr_h_stats[['Group', 'Obs_Name', 'Weight', 'Reduction_%', 'Final_Contribution_%']].to_markdown(index=False))
            arr_h_contrib = arr_h_stats['Final_Contribution_%'].sum()
            f.write(f"\n\nTotal ARR-H contribution to final phi: {arr_h_contrib:.2f}%\n\n")
        else:
            f.write("No ARR-H observations found.\n\n")

        f.write("## All Observation Groups\n\n")
        f.write(stats_df[['Type', 'Group', 'Obs_Name', 'Weight', 'Reduction_%', 'Final_Contribution_%']].to_markdown(index=False))

    print(f"Markdown summary saved to: {md_path}")

    print("\n" + "="*100)
    print("ANALYSIS COMPLETE!")
    print("="*100)

    return stats_df


if __name__ == '__main__':
    # Get model name from command line or use default
    if len(sys.argv) > 1:
        model_name = sys.argv[1]
    else:
        model_name = DEFAULT_MODEL_NAME
        print(f"No model name provided, using MODEL_NAME from setup.py: {model_name}")

    analyze_phi(model_name)
