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
import argparse

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

# Import setup to get model name
from setup import MODEL_NAME as DEFAULT_MODEL_NAME
from setup import REINFLATE_ITERS as DEFAULT_REINFLATE_LIST

def get_base_obgnme(obgnme):
    """Remove leading numbers from observation group names like '16obg4' -> 'obg4'"""
    match = re.match(r'^(\d+|greater_\d+|less_\d+)(.*)$', str(obgnme))
    if match and match.group(2):
        return match.group(2)
    return obgnme


def analyze_phi(run_name, model_name, reinflate_iters, last_iter_idx=None, selected_iter=None):
    """
    Analyze phi results for a given model run.

    Args:
        model_name: Name of the model run (e.g., 'local_run33')
    """
    print("="*100)
    print(f"PEST++ IES PHI ANALYSIS: {run_name} - {model_name}")
    print("="*100)

    print(os.getcwd())

    # Set paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_phi_analysis')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

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
            # Calculate reduction as percentage of final phi to allow for >100% reductions
            reduction = ((initial - final) / final * 100) if final > 0 else np.nan
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

    width = 12  # inches
    height = 0.6 * width  # Maintain aspect ratio
    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(
        3, 3, 
        hspace=0.4, wspace=0.2
        )

    # Plot A: Normalized phi evolution with line chart and percentile shading
    # Mark reinflation iterations first (for later use)

    ax1 = fig.add_subplot(gs[0, :])

    # Normalize all realizations to initial phi (iteration 0)
    # Use realization IDs to match across iterations
    iterations = sorted(phi_group['iteration'].unique())

    # Get initial phi for each realization (from iteration 0)
    init_data = phi_group[phi_group['iteration'] == 0].copy()
    init_data['total_phi'] = init_data[obs_group_cols].sum(axis=1)
    init_phi_map = dict(zip(init_data['par_realization'], init_data['total_phi']))

    # Calculate percentiles for each iteration
    percentile_data = {iter_num: [] for iter_num in iterations}

    for iter_num in iterations:
        iter_data = phi_group[phi_group['iteration'] == iter_num].copy()
        iter_data['total_phi'] = iter_data[obs_group_cols].sum(axis=1)

        # Normalize each realization by its initial phi
        for idx, row in iter_data.iterrows():
            real_id = row['par_realization']
            if real_id in init_phi_map and init_phi_map[real_id] > 0:
                normalized = row['total_phi'] / init_phi_map[real_id]
                percentile_data[iter_num].append(normalized)

    # Calculate median, IQR, and 5-95th percentile
    median_vals = [np.median(percentile_data[i]) if len(percentile_data[i]) > 0 else np.nan for i in iterations]
    q25_vals = [np.percentile(percentile_data[i], 25) if len(percentile_data[i]) > 0 else np.nan for i in iterations]
    q75_vals = [np.percentile(percentile_data[i], 75) if len(percentile_data[i]) > 0 else np.nan for i in iterations]
    p05_vals = [np.percentile(percentile_data[i], 5) if len(percentile_data[i]) > 0 else np.nan for i in iterations]
    p95_vals = [np.percentile(percentile_data[i], 95) if len(percentile_data[i]) > 0 else np.nan for i in iterations]

    # Plot shaded regions
    ax1.fill_between(iterations, p05_vals, p95_vals, alpha=0.2, color='blue', label='5-95th percentile')
    ax1.fill_between(iterations, q25_vals, q75_vals, alpha=0.3, color='blue', label='IQR (25-75th)')

    # Plot median line
    ax1.plot(iterations, median_vals, 'o-', linewidth=2, markersize=6, color='darkblue',
             label='Median', zorder=3)

    # Mark reinflation iterations
    for iter_num in reinflate_iters:
        ax1.axvline(x=iter_num, color='red', linestyle='--', linewidth=1, alpha=0.5, zorder=1)
    if len(reinflate_iters) > 0:
        ax1.plot([], [], 'r--', label='Reinflation', linewidth=1, alpha=0.5)

    ax1.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax1.set_xlabel('Iteration', fontsize=8)
    ax1.set_ylabel('Normalized φ', fontsize=8)
    ax1.set_title('A. Normalized φ decrease', fontsize=8, loc='left', fontweight='bold')
    ax1.legend(fontsize=6, loc='lower left')
    ax1.grid(True, which='both', alpha=0.3, axis='both')
    ax1.set_yscale('log')
    ax1.set_xlim(0, last_iter_idx)
    ax1.tick_params(labelsize=8)
    ax1.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Define function to convert observation names to prettier labels
    def get_obs_label(obs_name):
        """Convert observation names to prettier labels"""
        obs_lower = str(obs_name).lower()
        # Flux groups
        if 'arr-awq' in obs_lower or 'awq' in obs_lower:
            return 'Awanui flux'
        elif 'arr-spq' in obs_lower or 'arr-sprq' in obs_lower or 'spq' in obs_lower:
            return 'Spring flux'
        elif 'arr-confq' in obs_lower or 'confq' in obs_lower:
            return 'Confined flux'
        elif 'arr-pwq' in obs_lower or 'pwq' in obs_lower:
            return 'Poukawa flux'
        # Heads
        elif 'arr-h' in obs_lower:
            return 'Model heads'
        elif 'ts-heads' in obs_lower:
            return 'Timeseries heads'
        # Recession
        elif 'recession' in obs_lower and 'budget' not in obs_lower:
            return 'Recession rate'
        # Budget
        elif 'budget' in obs_lower:
            if 'confined' in obs_lower:
                return 'Budget: Confined'
            elif 'inflow' in obs_lower:
                return 'Budget: Rch+Conf'
            elif 'sw' in obs_lower:
                return 'Budget: SW'
            else:
                return obs_name  # Keep original for others
        else:
            return obs_name

    # Define color scheme function with darker colors as requested
    def get_obs_color(obs_name):
        """Assign colors based on observation subcategory"""
        obs_lower = obs_name.lower()
        # Flux groups - shades of red/pink
        if 'arr-awq' in obs_lower or 'awq' in obs_lower:
            return '#e74c3c'  # Red for Awanui
        elif 'arr-spq' in obs_lower or 'arr-sprq' in obs_lower or 'spq' in obs_lower:
            return '#8B0000'  # Darker red for springs (DarkRed)
        elif 'arr-confq' in obs_lower or 'confq' in obs_lower:
            return '#e67e22'  # Orange for confined arr
        elif 'arr-pwq' in obs_lower or 'pwq' in obs_lower:
            return '#d35400'  # Dark orange for Poukawa
        # Heads groups - shades of blue/green
        elif 'arr-h' in obs_lower:
            return '#3498db'  # Blue for arr-h
        elif 'ts-heads' in obs_lower:
            return '#1a5490'  # Darker blue for ts-heads
        # Recession - green
        elif 'recession' in obs_lower:
            return '#27ae60'  # Darker green for recession
        # Budget - darker purple/brown
        elif 'budget' in obs_lower:
            if 'confined' in obs_lower:
                return '#6a1b9a'  # Darker purple for confined budget
            elif 'inflow' in obs_lower:
                return '#9b59b6'  # Medium purple for inflow
            elif 'sw' in obs_lower:
                return '#8e44ad'  # Purple for sw
            else:
                return '#7b68ee'  # Medium slate blue for other budget
        else:
            return '#95a5a6'  # Gray for other

    # Plots B-D: Phi evolution by grouped observation types (normalized)
    # Filter to show only specific observations as requested:
    # B) Fluxes: arr-awq, arr-pwq, arr-confq, arr-sprq
    # C) Heads: one ts-heads (pk4), arr-h, one recession (pk4)
    # D) Budget: budget-confined, budget-inflow, budget-sw

    flux_groups_new = {}
    heads_groups_new = {}
    budget_groups_new = {}

    # Track if we've already added one ts-heads and one recession:pk4
    ts_heads_added = False
    recession_pk4_added = False

    for group, oname in obs_type_map.items():
        oname_lower = str(oname).lower()
        # Flux groups
        if any(x in oname_lower for x in ['arr-awq', 'arr-pwq', 'arr-confq', 'arr-spq', 'arr-sprq']):
            flux_groups_new[group] = oname
        # Heads groups - filter to only one ts-heads and recession:pk4
        elif any(x in oname_lower for x in ['arr-h']):
            heads_groups_new[group] = oname
        elif 'ts-heads' in oname_lower and not ts_heads_added:
            if 'pk4' in oname_lower and 'diff' not in oname_lower:
                heads_groups_new[group] = oname
                ts_heads_added = True
        elif 'recession' in oname_lower and 'budget' not in oname_lower and not recession_pk4_added:
            if 'pk4' in oname_lower and 'diff' not in oname_lower:
                heads_groups_new[group] = oname
                recession_pk4_added = True
        # Budget groups (budget-confined, budget-inflow, budget-sw)
        elif 'budget' in oname_lower:
            if any(x in oname_lower for x in ['confined', 'inflow', 'sw']):
                budget_groups_new[group] = oname

    group_configs = [
        (flux_groups_new, 'B. Fluxes', 'o', gs[1, 0]),
        (heads_groups_new, 'C. Heads', 's', gs[1, 1]),
        (budget_groups_new, 'D. Budget', '^', gs[1, 2])
    ]

    # # Get last iteration for x-axis limits
    # last_iteration = int(phi_by_iter.index[-1])

    for idx, (groups, title, marker, subplot_loc) in enumerate(group_configs):
        ax = fig.add_subplot(subplot_loc)
        if len(groups) > 0:
            for group, oname in groups.items():
                if group in phi_by_iter.columns:
                    # Only plot groups with non-zero weights
                    weight = phi_factors[phi_factors['group'] == group]['weight'].values
                    weight_val = weight[0] if len(weight) > 0 else 0.0
                    if weight_val > 0:
                        # Normalize to initial value
                        initial = phi_by_iter.loc[0, group]
                        normalized = phi_by_iter[group] / initial if initial > 0 else phi_by_iter[group]
                        color = get_obs_color(oname)
                        label = get_obs_label(oname)
                        ax.plot(phi_by_iter.index, normalized, marker=marker,
                               linewidth=1.5, label=label, markersize=3, alpha=0.8, color=color)

            for iter_num in reinflate_iters:
                ax.axvline(x=iter_num, color='red', linestyle='--', linewidth=1, alpha=0.3)

            ax.legend(fontsize=6, loc='upper right')
        else:
            ax.text(0.5, 0.5, 'No observations', ha='center', va='center',
                   transform=ax.transAxes, fontsize=8, color='gray')

        ax.set_xlabel('Iteration', fontsize=8)
        if idx == 0:  # B. Fluxes
            ax.set_ylabel('Normalized φ', fontsize=8)
        else:
            # no label
            ax.set_ylabel('')
        ax.set_title(title, fontsize=8, loc='left', fontweight='bold')
        ax.set_yscale('log')

        # Set x-axis to show integers and limits
        ax.set_xlim(0, last_iter_idx)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

        # Grid settings
        ax.grid(True, which='major', alpha=0.3, axis='both')
        ax.tick_params(labelsize=8)

        # Add minor gridlines and tick marks for Plot C (Heads) on y-axis
        if idx == 1:  # C. Heads
            from matplotlib.ticker import LogLocator
            ax.yaxis.set_minor_locator(LogLocator(subs='auto'))
            ax.grid(True, which='minor', alpha=0.15, axis='y', linestyle=':')
            # Enable minor tick marks on y-axis
            ax.tick_params(axis='y', which='minor', length=3, width=0.5)

    # Bottom row: Two horizontal bar charts (plots E & F) sharing y-axis
    # Calculate statistics for second-to-last iteration (last_iter - 1)
    if last_iter_idx is None:
        last_iter_idx = phi_by_iter.index[-1]

    if selected_iter is None:
        selected_iter = last_iter_idx
    
    # second_last_iter_idx = phi_by_iter.index[-2] if len(phi_by_iter.index) > 1 else last_iter_idx

    # Recalculate stats for second-to-last iteration
    # Filter: only include one recession rate and one TS heads
    stats_second_last = []
    ts_heads_added = False
    recession_added = False

    for group in obs_group_cols:
        if group not in phi_by_iter.columns:
            continue

        initial = phi_by_iter.loc[0, group]
        second_last = phi_by_iter.loc[selected_iter, group]
        # Calculate reduction as percentage of final phi to allow for >100% reductions
        # E.g., (10^11 - 10^5) / 10^5 * 100 = very large percentage
        reduction = ((initial - second_last) / second_last * 100) if second_last > 0 else np.nan
        weight = phi_factors[phi_factors['group'] == group]['weight'].values
        weight_val = weight[0] if len(weight) > 0 else 0.0

        # Get observation name
        obs_name = obs_type_map.get(group, group)
        obs_lower = str(obs_name).lower()

        if weight_val > 0:  # Only include non-zero weights
            # Filter: only add one TS heads and one recession
            if 'ts-heads' in obs_lower:
                if ts_heads_added:
                    continue
                # Only add pk4 without diff
                if 'pk4' in obs_lower and 'diff' not in obs_lower:
                    ts_heads_added = True
                else:
                    continue
            elif 'recession' in obs_lower and 'budget' not in obs_lower:
                if recession_added:
                    continue
                # Only add pk4 without diff
                if 'pk4' in obs_lower and 'diff' not in obs_lower:
                    recession_added = True
                else:
                    continue

            stats_second_last.append({
                'Group': group,
                'Obs_Name': obs_name,
                'Weight': weight_val,
                'Final_Phi': second_last,
                'Reduction_%': reduction
            })

    stats_second_last_df = pd.DataFrame(stats_second_last)

    # Plot E: Final phi contributions by group (horizontal bar chart) - center position
    ax_e = fig.add_subplot(gs[2, 0:2])

    if len(stats_second_last_df) > 0:
        final_phi_df = stats_second_last_df.sort_values('Final_Phi', ascending=True)  # Ascending for horizontal
        # Apply prettier labels
        final_phi_df['Pretty_Name'] = final_phi_df['Obs_Name'].apply(get_obs_label)
        bar_colors = [get_obs_color(name) for name in final_phi_df['Obs_Name']]

        # Set grid behind bars
        ax_e.set_axisbelow(True)
        ax_e.grid(True, which='both', alpha=0.3, axis='x', zorder=0)

        bars = ax_e.barh(range(len(final_phi_df)), final_phi_df['Final_Phi'],
                       color=bar_colors, alpha=0.8, edgecolor='black', zorder=3)
        ax_e.set_yticks(range(len(final_phi_df)))
        ax_e.set_yticklabels(final_phi_df['Pretty_Name'], fontsize=8)
        ax_e.set_xlabel('φ', fontsize=8)
        ax_e.set_title(f'E. Final φ contributions by group (Iteration {selected_iter})',
                     fontsize=8, loc='left', fontweight='bold')
        ax_e.set_xscale('log')
        ax_e.tick_params(labelsize=8)
        # Ensure y-axis labels are visible
        plt.setp(ax_e.get_yticklabels(), visible=True)
    else:
        ax_e.text(0.5, 0.5, 'No non-zero weighted observations', ha='center', va='center',
                transform=ax_e.transAxes, fontsize=8, color='gray')

    # Plot F: Phi reduction from initial (horizontal bar chart) - right position, sharing y-axis
    ax_f = fig.add_subplot(gs[2, 2], sharey=ax_e)

    if len(stats_second_last_df) > 0:
        reduction_df = stats_second_last_df.sort_values('Final_Phi', ascending=True)  # Match E's order
        # Apply prettier labels
        reduction_df['Pretty_Name'] = reduction_df['Obs_Name'].apply(get_obs_label)
        bar_colors = [get_obs_color(name) for name in reduction_df['Obs_Name']]

        # Set grid behind bars
        ax_f.set_axisbelow(True)
        ax_f.grid(True, which='both', alpha=0.3, axis='x', zorder=0)

        bars = ax_f.barh(range(len(reduction_df)), reduction_df['Reduction_%'],
                       color=bar_colors, alpha=0.8, edgecolor='black', zorder=3)
        ax_f.set_yticks(range(len(reduction_df)))
        # Hide y-tick labels since we're sharing with ax_e
        # ax_f.set_yticklabels([])
        ax_f.set_xlabel('Reduction (%)', fontsize=8)
        ax_f.set_title(f'F. φ reduction from prior (Iteration {selected_iter})',
                     fontsize=8, loc='left', fontweight='bold')
        ax_f.set_xscale('log')
        ax_f.tick_params(labelsize=8)
        # Explicitly hide y-axis labels on F
        plt.setp(ax_f.get_yticklabels(), visible=False)
    else:
        ax_f.text(0.5, 0.5, 'No non-zero weighted observations', ha='center', va='center',
                transform=ax_f.transAxes, fontsize=8, color='gray')

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


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze PEST++ IES phi results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default=DEFAULT_MODEL_NAME,
        help=f'Model/run name for file prefixes (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )

    parser.add_argument(
        '--run-name', '-r',
        type=str,
        default=None,
        help='Run name for directory path (default: same as model_name)'
    )

    parser.add_argument(
        '--last-iter', '-l',
        type=int,
        default=None,
        help='Last iteration index for x-axis limit (default: auto-detect)'
    )

    parser.add_argument(
        '--selected-iter', '-s',
        type=int,
        default=None,
        help='Selected iteration for bar charts E & F (default: same as last-iter)'
    )

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    reinflate_iters = DEFAULT_REINFLATE_LIST

    r_iters = []
    start = 0
    for i in reinflate_iters:
        r_iters.append(start+i+1)
        start = start+i+1

    analyze_phi(run_name, model_name, r_iters, last_iter_idx=args.last_iter, selected_iter=args.selected_iter)
