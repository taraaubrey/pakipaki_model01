"""
Analyze parameter ensemble mean evolution across reinflation iterations.
Compares parameter means before/after reinflation events to understand
how ensemble re-inflation affects parameter estimates.

Usage:
    python scripts/i_analyze_param_reinflation.py [model_name]

    If model_name is not provided, uses MODEL_NAME from setup.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os
import sys
import glob
import pyemu

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

# Import setup and parameter categorization from h_analyze_params
from setup import MODEL_NAME as DEFAULT_MODEL_NAME
from h_analyze_params import categorize_parameter


def identify_reinflation_iters(phi_actual_file):
    """
    Identify iterations where reinflation occurred (phi increases).

    Args:
        phi_actual_file: Path to phi.actual.csv file

    Returns:
        list: Iteration numbers where reinflation occurred
    """
    phi = pd.read_csv(phi_actual_file)

    # Calculate phi changes
    phi['phi_change'] = phi['mean'].diff()

    # Reinflation occurs when phi increases significantly
    # Use a threshold: phi increases by >10x from previous iteration
    reinflate_iters = []
    for i in range(1, len(phi)):
        if phi.loc[i, 'mean'] > phi.loc[i-1, 'mean'] * 2:
            reinflate_iters.append(int(phi.loc[i, 'iteration']))

    return reinflate_iters


def load_parameter_means(par_file, pst):
    """
    Load parameter ensemble and calculate means by group.

    Args:
        par_file: Path to .par.csv file
        pst: pyemu.Pst object

    Returns:
        dict: {pargp: {'pname': str, 'mean': float, 'std': float, 'cv': float}}
    """
    # Load parameter ensemble
    par_df = pd.read_csv(par_file, index_col=0, low_memory=False)

    # Get parameter groups from PST
    param_data = pst.parameter_data.copy()

    # Calculate statistics by parameter group
    group_stats = {}

    for pargp in param_data['pargp'].unique():
        # Get all parameters in this group
        group_params = param_data[param_data['pargp'] == pargp].index.tolist()

        # Filter to parameters that exist in par_df
        group_params = [p for p in group_params if p in par_df.columns]

        if len(group_params) == 0:
            continue

        # Get values for this group across all realizations
        group_values = par_df[group_params].values.flatten()

        # Calculate statistics
        mean_val = np.mean(group_values)
        std_val = np.std(group_values)
        cv_val = (std_val / mean_val * 100) if mean_val != 0 else np.nan

        # Get pname from first parameter
        if 'pname' in param_data.columns:
            pname = param_data.loc[group_params[0], 'pname']
        else:
            # Extract from parameter name
            first_param = group_params[0]
            if 'pname:' in first_param:
                pname = first_param.split('pname:')[1].split('_')[0]
            else:
                pname = pargp

        group_stats[pargp] = {
            'pname': pname,
            'mean': mean_val,
            'std': std_val,
            'cv': cv_val,
            'n_params': len(group_params)
        }

    return group_stats


def analyze_reinflation_impact(model_name):
    """
    Analyze parameter mean evolution across reinflation iterations.

    Args:
        model_name: Name of the model run (e.g., 'local_run32')
    """
    print("="*100)
    print(f"PEST++ IES PARAMETER REINFLATION ANALYSIS: {model_name}")
    print("="*100)

    # Set paths
    master_dir = os.path.join('models', model_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', model_name)

    if not os.path.exists(master_dir):
        print(f"ERROR: Directory not found: {master_dir}")
        return

    # Load PST file
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    if not os.path.exists(pst_file):
        print(f"ERROR: PST file not found: {pst_file}")
        return

    pst = pyemu.Pst(pst_file)
    print(f"\nLoaded PST file: {pst_file}")
    print(f"Total parameters: {pst.npar}")
    print(f"Parameter groups: {len(pst.par_groups)}")

    # Identify reinflation iterations
    phi_file = os.path.join(master_dir, f'{model_name}.phi.actual.csv')
    if not os.path.exists(phi_file):
        print(f"ERROR: Phi file not found: {phi_file}")
        return

    reinflate_iters = identify_reinflation_iters(phi_file)
    print(f"\nReinflation iterations detected: {reinflate_iters}")

    # Find all parameter files
    par_files = sorted(glob.glob(os.path.join(master_dir, f'{model_name}.*.par.csv')))

    # Extract iteration numbers
    iterations = []
    for pf in par_files:
        basename = os.path.basename(pf)
        iter_num = int(basename.replace('.par.csv', '').split('.')[-1])
        iterations.append(iter_num)

    iterations = sorted(iterations)
    print(f"\nParameter files found for iterations: {iterations}")

    # Define key iterations to analyze
    # Include: iter 0, iterations before/after each reinflation, and final
    key_iters = [0]  # Always include initial
    for r_iter in reinflate_iters:
        if r_iter - 1 > 0:
            key_iters.append(r_iter - 1)  # Before reinflation
        key_iters.append(r_iter)  # Reinflation
        if r_iter + 1 in iterations:
            key_iters.append(r_iter + 1)  # After reinflation

    if iterations[-1] not in key_iters:
        key_iters.append(iterations[-1])  # Final iteration

    key_iters = sorted(set(key_iters))
    print(f"\nAnalyzing key iterations: {key_iters}")

    # Load parameter means for all iterations
    print("\nLoading parameter ensemble data...")
    all_stats = {}

    for iter_num in iterations:
        par_file = os.path.join(master_dir, f'{model_name}.{iter_num}.par.csv')
        if os.path.exists(par_file):
            all_stats[iter_num] = load_parameter_means(par_file, pst)
            print(f"  Iteration {iter_num}: {len(all_stats[iter_num])} parameter groups")

    # Organize data by parameter group
    print("\nOrganizing data by parameter group...")
    param_groups = list(all_stats[0].keys())

    # Create evolution dataframe for each group
    evolution_data = []
    for pargp in param_groups:
        pname = all_stats[0][pargp]['pname']
        category, subcategory = categorize_parameter(pname)

        for iter_num in iterations:
            if pargp in all_stats[iter_num]:
                evolution_data.append({
                    'Iteration': iter_num,
                    'Group': pargp,
                    'Pname': pname,
                    'Category': category,
                    'Subcategory': subcategory,
                    'Mean': all_stats[iter_num][pargp]['mean'],
                    'Std': all_stats[iter_num][pargp]['std'],
                    'CV': all_stats[iter_num][pargp]['cv'],
                    'N_Params': all_stats[iter_num][pargp]['n_params'],
                    'Is_Reinflation': iter_num in reinflate_iters
                })

    evolution_df = pd.DataFrame(evolution_data)

    # Calculate mean change at reinflation points
    print("\n" + "="*100)
    print("PARAMETER MEAN CHANGES AT REINFLATION ITERATIONS")
    print("="*100)

    for r_iter in reinflate_iters:
        print(f"\nReinflation at iteration {r_iter}:")
        print(f"{'Parameter':<30s} {'Before':>15s} {'After':>15s} {'Change %':>15s}")
        print("-"*80)

        before_iter = r_iter - 1
        after_iter = r_iter

        for pargp in param_groups:
            if pargp in all_stats[before_iter] and pargp in all_stats[after_iter]:
                pname = all_stats[before_iter][pargp]['pname']
                mean_before = all_stats[before_iter][pargp]['mean']
                mean_after = all_stats[after_iter][pargp]['mean']
                change_pct = ((mean_after - mean_before) / mean_before * 100) if mean_before != 0 else np.nan

                print(f"{pname:<30s} {mean_before:>15.6f} {mean_after:>15.6f} {change_pct:>14.2f}%")

    # Load phi data for comparison
    phi_actual = pd.read_csv(phi_file)

    # Load PCS data if available
    pcs_files = sorted(glob.glob(os.path.join(master_dir, f'{model_name}.*.pcs.csv')))
    pcs_data = []
    if len(pcs_files) > 0:
        for pcs_file in pcs_files:
            basename = os.path.basename(pcs_file)
            # Skip files that don't match the iteration pattern
            try:
                iter_num = int(basename.replace('.pcs.csv', '').split('.')[-1])
            except ValueError:
                # Skip files like .reinflate.pcs.csv
                continue
            pcs = pd.read_csv(pcs_file)
            pcs['iteration'] = iter_num
            pcs_data.append(pcs)
        if len(pcs_data) > 0:
            pcs_df = pd.concat(pcs_data, ignore_index=True)
        else:
            pcs_df = None
    else:
        pcs_df = None

    # Create visualizations
    print("\n" + "="*100)
    print("CREATING VISUALIZATIONS...")
    print("="*100)

    fig = plt.figure(figsize=(20, 16))
    gs = fig.add_gridspec(5, 3, hspace=0.35, wspace=0.3)

    # Color scheme for parameter categories
    category_colors = {
        'Hydraulic Conductivity': '#e74c3c',
        'Storage': '#3498db',
        'Recharge': '#2ecc71',
        'GHB Conductance': '#f39c12',
        'GHB Head': '#9b59b6',
        'GHB Other': '#95a5a6',
        'Other': '#34495e'
    }

    # Plot 1: Parameter means evolution across all iterations (by category)
    # GHB Conductance on left log scale, others on right log scale
    ax1 = fig.add_subplot(gs[0, :])

    # Separate GHB Conductance from other categories
    ghb_cond_categories = ['GHB Conductance']
    other_categories = [cat for cat in evolution_df['Category'].unique()
                       if cat not in ghb_cond_categories]

    # Plot GHB Conductance on primary (log) axis
    for category in ghb_cond_categories:
        if category in evolution_df['Category'].values:
            cat_data = evolution_df[evolution_df['Category'] == category]
            cat_mean = cat_data.groupby('Iteration')['Mean'].mean()
            color = category_colors.get(category, '#95a5a6')
            ax1.plot(cat_mean.index, cat_mean.values, 'o-', linewidth=2.5,
                    markersize=6, label=category, color=color, alpha=0.8)

    ax1.set_yscale('log')
    ax1.set_ylabel('GHB Conductance Mean (log)', fontsize=11, fontweight='bold', color='#f39c12')
    ax1.set_xlabel('Iteration', fontsize=12, fontweight='bold')
    ax1.tick_params(axis='y', labelcolor='#f39c12')

    # Create secondary axis for other categories (also log scale)
    ax1b = ax1.twinx()
    for category in other_categories:
        cat_data = evolution_df[evolution_df['Category'] == category]
        cat_mean = cat_data.groupby('Iteration')['Mean'].mean()
        color = category_colors.get(category, '#95a5a6')
        ax1b.plot(cat_mean.index, cat_mean.values, 's-', linewidth=2,
                 markersize=5, label=category, color=color, alpha=0.7)

    ax1b.set_yscale('log')
    ax1b.set_ylabel('Other Parameter Means (log)', fontsize=11, fontweight='bold')

    # Mark reinflation iterations on both axes
    for r_iter in reinflate_iters:
        ax1.axvline(x=r_iter, color='red', linestyle='--', linewidth=2, alpha=0.5)

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax1b.get_legend_handles_labels()
    if len(reinflate_iters) > 0:
        lines1.append(plt.Line2D([0], [0], color='red', linestyle='--', linewidth=2, alpha=0.5))
        labels1.append('Reinflation')
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=9, ncol=2)

    ax1.set_title(f'Parameter Mean Evolution Across Reinflation: {model_name}',
                 fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Plot 2: CV evolution by category
    ax2 = fig.add_subplot(gs[1, 0])
    for category in evolution_df['Category'].unique():
        cat_data = evolution_df[evolution_df['Category'] == category]
        cat_cv = cat_data.groupby('Iteration')['CV'].mean()
        color = category_colors.get(category, '#95a5a6')
        ax2.plot(cat_cv.index, cat_cv.values, 'o-', linewidth=2,
                markersize=5, label=category, color=color, alpha=0.8)

    for r_iter in reinflate_iters:
        ax2.axvline(x=r_iter, color='red', linestyle='--', linewidth=1.5, alpha=0.5)

    ax2.set_xlabel('Iteration', fontsize=11)
    ax2.set_ylabel('Coefficient of Variation (%)', fontsize=11)
    ax2.set_title('Ensemble Spread (CV) Evolution', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=8, loc='best')
    ax2.grid(True, alpha=0.3)

    # Plot 3: Phi evolution
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.plot(phi_actual['iteration'], phi_actual['mean'], 'o-', linewidth=2.5,
            markersize=7, color='black', label='Mean Phi', zorder=3)
    ax3.fill_between(phi_actual['iteration'],
                     phi_actual['mean'] - phi_actual['standard_deviation'],
                     phi_actual['mean'] + phi_actual['standard_deviation'],
                     alpha=0.2, color='black', label='±1 Std')

    for r_iter in reinflate_iters:
        ax3.axvline(x=r_iter, color='red', linestyle='--', linewidth=2, alpha=0.5)

    ax3.set_xlabel('Iteration', fontsize=11)
    ax3.set_ylabel('Total Phi', fontsize=11)
    ax3.set_title('Phi Evolution', fontsize=12, fontweight='bold')
    ax3.set_yscale('log')
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.3)

    # Plot 4: Mean change magnitude at reinflation vs CV
    if pcs_df is not None:
        ax4 = fig.add_subplot(gs[1, 2])

        # Calculate CV change relative to each reinflation period
        # Determine period start iterations (0, first_reinflate, second_reinflate, etc.)
        period_starts = [0] + sorted(reinflate_iters)

        all_mean_changes = []
        all_cv_changes = []
        all_iterations = []
        all_categories = []
        all_periods = []  # Track which period each point belongs to

        for iter_num in iterations:
            # Find which reinflation period this iteration belongs to
            period_start = 0
            period_idx = 0
            for idx, ps in enumerate(reversed(period_starts)):
                if iter_num >= ps:
                    period_start = ps
                    period_idx = len(period_starts) - 1 - idx
                    break

            if iter_num in pcs_df['iteration'].values:
                pcs_iter = pcs_df[pcs_df['iteration'] == iter_num]

                for _, row in pcs_iter.iterrows():
                    pargp = row['group']
                    if pargp in evolution_df['Group'].values:
                        mean_change = row['mean_change']

                        # Get CV at current iteration and at period start
                        current_cv = evolution_df[(evolution_df['Group'] == pargp) &
                                                 (evolution_df['Iteration'] == iter_num)]['CV'].values
                        period_start_cv = evolution_df[(evolution_df['Group'] == pargp) &
                                                      (evolution_df['Iteration'] == period_start)]['CV'].values

                        if len(current_cv) > 0 and len(period_start_cv) > 0:
                            cv_change = current_cv[0] - period_start_cv[0]
                        else:
                            cv_change = 0

                        # Get category for this group
                        category = evolution_df[evolution_df['Group'] == pargp]['Category'].iloc[0]

                        all_mean_changes.append(mean_change)
                        all_cv_changes.append(cv_change)
                        all_iterations.append(iter_num)
                        all_categories.append(category)
                        all_periods.append(period_idx)

        # Define colors for each period
        period_colors = ['#3498db', '#e74c3c', '#2ecc71', '#f39c12', '#9b59b6', '#1abc9c']
        period_labels = []
        for i in range(len(period_starts)):
            if i < len(period_starts) - 1:
                if i == 0:
                    period_labels.append(f'Period {i+1} (iter 0-{period_starts[i+1]-1})')
                else:
                    next_start = period_starts[i+1] if i+1 < len(period_starts) else iterations[-1]+1
                    period_labels.append(f'Period {i+1} (iter {period_starts[i]}-{next_start-1})')
            else:
                period_labels.append(f'Period {i+1} (iter {period_starts[i]}-{iterations[-1]})')

        # Define marker styles for each category
        category_markers = {
            'Hydraulic Conductivity': 'o',
            'Storage': 's',
            'Recharge': '^',
            'GHB Conductance': 'D',
            'GHB Head': 'v',
            'GHB Other': 'p',
            'Other': 'h'
        }

        # Plot scatter points colored by period with markers by category
        unique_categories = list(set(all_categories))
        unique_periods = list(set(all_periods))

        for period_idx in unique_periods:
            period_mask = [p == period_idx for p in all_periods]
            period_mean = [all_mean_changes[i] for i in range(len(all_mean_changes)) if period_mask[i]]
            period_cv = [all_cv_changes[i] for i in range(len(all_cv_changes)) if period_mask[i]]
            period_cats = [all_categories[i] for i in range(len(all_categories)) if period_mask[i]]

            color = period_colors[period_idx % len(period_colors)]

            # Plot points grouped by category to get correct markers
            for category in unique_categories:
                cat_mask = [c == category for c in period_cats]
                cat_mean = [period_mean[i] for i in range(len(period_mean)) if cat_mask[i]]
                cat_cv = [period_cv[i] for i in range(len(period_cv)) if cat_mask[i]]

                if len(cat_mean) > 0:
                    marker = category_markers.get(category, 'o')
                    ax4.scatter(cat_mean, cat_cv,
                              c=color, marker=marker, s=30, alpha=0.6,
                              edgecolors='black', linewidth=0.5, zorder=2)

        # Add category markers to legend
        from matplotlib.lines import Line2D
        legend_elements = []

        # Add period colors
        for i, label in enumerate(period_labels):
            legend_elements.append(Line2D([0], [0], marker='o', color='w',
                                        markerfacecolor=period_colors[i % len(period_colors)],
                                        markersize=8, label=label,
                                        markeredgecolor='black', markeredgewidth=0.5))

        # Add separator
        legend_elements.append(Line2D([0], [0], linestyle='none', label=''))

        # Add category markers
        for category, marker in category_markers.items():
            if category in unique_categories:
                legend_elements.append(Line2D([0], [0], marker=marker, color='w',
                                            markerfacecolor='gray', markersize=7,
                                            label=category, markeredgecolor='black',
                                            markeredgewidth=0.5))

        # Mark reinflation iterations with special markers
        for r_iter in reinflate_iters:
            if r_iter in all_iterations:
                r_indices = [i for i, it in enumerate(all_iterations) if it == r_iter]
                r_mean_changes = [all_mean_changes[i] for i in r_indices]
                r_cv_changes = [all_cv_changes[i] for i in r_indices]
                ax4.scatter(r_mean_changes, r_cv_changes,
                          marker='X', s=75, c='red', alpha=0.9,
                          edgecolors='black', linewidth=1.5, zorder=10)
                legend_elements.append(Line2D([0], [0], marker='X', color='w',
                                            markerfacecolor='red', markersize=9,
                                            label=f'Iter {r_iter} (reinflate)',
                                            markeredgecolor='black', markeredgewidth=1.5))

        # Set log scales for both axes
        ax4.set_xscale('log')
        ax4.set_yscale('symlog')  # Use symlog for y-axis to handle negative values

        ax4.axhline(y=0, color='gray', linestyle=':', linewidth=1)
        ax4.set_xlabel('Mean Parameter Change (log scale)', fontsize=11)
        ax4.set_ylabel('CV Change from Period Start (symlog)', fontsize=11)
        ax4.set_title('Parameter Change vs CV Change (color=period, marker=category)',
                     fontsize=11, fontweight='bold')
        ax4.legend(handles=legend_elements, fontsize=6, loc='best', ncol=2)
        ax4.grid(True, alpha=0.3, which='both')

    # Plots 5-7: Individual parameter group evolution for selected groups
    # Select top 3 categories with most parameters
    top_categories = evolution_df.groupby('Category')['N_Params'].sum().nlargest(3).index.tolist()

    # Define marker styles for different parameterization types
    marker_map = {
        'hg': 'o',      # head group (time-varying)
        'cn': 's',      # constant
        'gr': '^',      # grid
        'pp': 'D',      # pilot points
        'other': 'v'
    }

    # Define color mapping for K parameterization types
    k_color_map = {
        'K-grid': '#e74c3c',      # Red
        'K-pilot': '#3498db',     # Blue
        'K-constant': '#2ecc71'   # Green
    }

    for idx, category in enumerate(top_categories):
        ax = fig.add_subplot(gs[2, idx])
        cat_data = evolution_df[evolution_df['Category'] == category]

        # For Hydraulic Conductivity, use specific colors for grid/pilot/constant
        if category == 'Hydraulic Conductivity':
            for subcategory in cat_data['Subcategory'].unique():
                subcat_data = cat_data[cat_data['Subcategory'] == subcategory]

                # Get color and marker based on subcategory
                color = k_color_map.get(subcategory, '#95a5a6')
                # Extract suffix for marker
                parts = subcategory.rsplit('-', 1)
                suffix = parts[1] if len(parts) == 2 else 'other'
                marker = marker_map.get(suffix, marker_map['other'])

                ax.plot(subcat_data['Iteration'], subcat_data['Mean'],
                       marker=marker, linestyle='-', linewidth=2, markersize=5,
                       label=subcategory, color=color, alpha=0.8)
        else:
            # For other categories, group by base parameter
            base_groups = {}
            for subcategory in cat_data['Subcategory'].unique():
                # Extract base name (everything before last hyphen)
                parts = subcategory.rsplit('-', 1)
                if len(parts) == 2:
                    base_name = parts[0]
                    suffix = parts[1]
                else:
                    base_name = subcategory
                    suffix = 'other'

                if base_name not in base_groups:
                    base_groups[base_name] = []
                base_groups[base_name].append((subcategory, suffix))

            # Assign colors to base groups
            base_names = list(base_groups.keys())
            colors = plt.cm.tab10(np.linspace(0, 1, len(base_names)))

            # Plot each subcategory with color by base group, marker by suffix
            for base_idx, (base_name, subcat_list) in enumerate(base_groups.items()):
                color = colors[base_idx]

                for subcategory, suffix in subcat_list:
                    subcat_data = cat_data[cat_data['Subcategory'] == subcategory]
                    marker = marker_map.get(suffix, marker_map['other'])

                    ax.plot(subcat_data['Iteration'], subcat_data['Mean'],
                           marker=marker, linestyle='-', linewidth=2, markersize=5,
                           label=subcategory, color=color, alpha=0.8)

        for r_iter in reinflate_iters:
            ax.axvline(x=r_iter, color='red', linestyle='--', linewidth=1, alpha=0.5)

        ax.set_xlabel('Iteration', fontsize=10)
        ax.set_ylabel('Parameter Mean', fontsize=10)
        ax.set_title(f'{category} Evolution', fontsize=11, fontweight='bold')
        ax.legend(fontsize=6, loc='best', ncol=1)
        ax.grid(True, alpha=0.3)

    # Plot 8: CV at key iterations (before/after reinflation)
    ax8 = fig.add_subplot(gs[3, :2])

    key_iter_data = evolution_df[evolution_df['Iteration'].isin(key_iters)]

    # Create grouped bar chart
    categories = key_iter_data['Category'].unique()
    x = np.arange(len(categories))
    width = 0.8 / len(key_iters)

    for i, iter_num in enumerate(key_iters):
        iter_data = key_iter_data[key_iter_data['Iteration'] == iter_num]
        cv_by_cat = [iter_data[iter_data['Category'] == cat]['CV'].mean()
                     for cat in categories]

        is_reinflate = iter_num in reinflate_iters
        color = 'red' if is_reinflate else 'blue'
        alpha = 0.9 if is_reinflate else 0.6

        ax8.bar(x + i * width, cv_by_cat, width, label=f'Iter {iter_num}',
               alpha=alpha, edgecolor='black', linewidth=0.5,
               color=color if is_reinflate else None)

    ax8.set_xlabel('Parameter Category', fontsize=11)
    ax8.set_ylabel('Coefficient of Variation (%)', fontsize=11)
    ax8.set_title('CV at Key Iterations (Reinflation Comparison)', fontsize=12, fontweight='bold')
    ax8.set_xticks(x + width * (len(key_iters) - 1) / 2)
    ax8.set_xticklabels(categories, rotation=45, ha='right', fontsize=9)
    ax8.legend(fontsize=8, loc='best', ncol=min(4, len(key_iters)))
    ax8.grid(True, alpha=0.3, axis='y')

    # Plot 9: Parameter mean percent change from initial
    ax9 = fig.add_subplot(gs[3, 2])

    # Calculate percent change from iteration 0
    for category in evolution_df['Category'].unique():
        cat_data = evolution_df[evolution_df['Category'] == category]

        # Average initial mean for this category
        initial_mean = cat_data[cat_data['Iteration'] == 0]['Mean'].mean()

        # Calculate percent change for each iteration
        pct_changes = []
        iter_nums = []
        for iter_num in iterations:
            iter_mean = cat_data[cat_data['Iteration'] == iter_num]['Mean'].mean()
            pct_change = ((iter_mean - initial_mean) / initial_mean * 100) if initial_mean != 0 else 0
            pct_changes.append(pct_change)
            iter_nums.append(iter_num)

        color = category_colors.get(category, '#95a5a6')
        ax9.plot(iter_nums, pct_changes, 'o-', linewidth=2, markersize=5,
                label=category, color=color, alpha=0.8)

    for r_iter in reinflate_iters:
        ax9.axvline(x=r_iter, color='red', linestyle='--', linewidth=1.5, alpha=0.5)

    ax9.axhline(y=0, color='gray', linestyle=':', linewidth=1)
    ax9.set_xlabel('Iteration', fontsize=11)
    ax9.set_ylabel('Change from Initial (%)', fontsize=11)
    ax9.set_title('Parameter Mean % Change from Initial', fontsize=12, fontweight='bold')
    ax9.legend(fontsize=8, loc='best')
    ax9.grid(True, alpha=0.3)

    # Plot 10: Summary statistics table
    ax10 = fig.add_subplot(gs[4, :])
    ax10.axis('off')

    # Create summary table
    summary_rows = []
    summary_rows.append(['Category', 'Initial Mean', 'Final Mean', 'Change %',
                        'Initial CV', 'Final CV', 'CV Change'])

    for category in evolution_df['Category'].unique():
        cat_data = evolution_df[evolution_df['Category'] == category]

        initial = cat_data[cat_data['Iteration'] == 0]
        final = cat_data[cat_data['Iteration'] == iterations[-1]]

        init_mean = initial['Mean'].mean()
        final_mean = final['Mean'].mean()
        mean_change = ((final_mean - init_mean) / init_mean * 100) if init_mean != 0 else 0

        init_cv = initial['CV'].mean()
        final_cv = final['CV'].mean()
        cv_change = final_cv - init_cv

        summary_rows.append([
            category,
            f'{init_mean:.4f}',
            f'{final_mean:.4f}',
            f'{mean_change:+.2f}%',
            f'{init_cv:.2f}',
            f'{final_cv:.2f}',
            f'{cv_change:+.2f}'
        ])

    table = ax10.table(cellText=summary_rows, cellLoc='center', loc='center',
                      bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)

    # Style header row
    for i in range(len(summary_rows[0])):
        table[(0, i)].set_facecolor('#3498db')
        table[(0, i)].set_text_props(weight='bold', color='white')

    # Alternate row colors
    for i in range(1, len(summary_rows)):
        color = '#ecf0f1' if i % 2 == 0 else 'white'
        for j in range(len(summary_rows[0])):
            table[(i, j)].set_facecolor(color)

    # Save figure
    fig_path = os.path.join(output_dir, f'{model_name}_reinflation_analysis.png')
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"Figure saved to: {fig_path}")
    plt.close()

    # Save evolution data to CSV
    csv_path = os.path.join(output_dir, f'{model_name}_param_evolution.csv')
    evolution_df.to_csv(csv_path, index=False)
    print(f"Evolution data saved to: {csv_path}")

    # Create markdown summary
    md_path = os.path.join(output_dir, f'{model_name}_REINFLATION_ANALYSIS.md')
    with open(md_path, 'w') as f:
        f.write(f"# Parameter Reinflation Analysis: {model_name}\n\n")

        f.write("## Overview\n\n")
        f.write(f"Total iterations: {len(iterations)}\n\n")
        f.write(f"Reinflation iterations: {', '.join(map(str, reinflate_iters))}\n\n")
        f.write(f"Parameter groups analyzed: {len(param_groups)}\n\n")

        f.write("## Parameter Mean Changes at Reinflation\n\n")
        for r_iter in reinflate_iters:
            f.write(f"### Iteration {r_iter} (Reinflation)\n\n")

            before_iter = r_iter - 1
            change_data = []

            for pargp in param_groups:
                if pargp in all_stats[before_iter] and pargp in all_stats[r_iter]:
                    pname = all_stats[before_iter][pargp]['pname']
                    mean_before = all_stats[before_iter][pargp]['mean']
                    mean_after = all_stats[r_iter][pargp]['mean']
                    cv_before = all_stats[before_iter][pargp]['cv']
                    cv_after = all_stats[r_iter][pargp]['cv']
                    mean_change = ((mean_after - mean_before) / mean_before * 100) if mean_before != 0 else np.nan
                    cv_change = cv_after - cv_before

                    change_data.append({
                        'Parameter': pname,
                        'Mean_Before': mean_before,
                        'Mean_After': mean_after,
                        'Mean_Change_%': mean_change,
                        'CV_Before': cv_before,
                        'CV_After': cv_after,
                        'CV_Change': cv_change
                    })

            change_df = pd.DataFrame(change_data)
            change_df = change_df.sort_values('Mean_Change_%', key=abs, ascending=False)
            f.write(change_df.to_markdown(index=False))
            f.write("\n\n")

        f.write("## Summary Statistics\n\n")
        f.write(pd.DataFrame(summary_rows[1:], columns=summary_rows[0]).to_markdown(index=False))
        f.write("\n\n")

        f.write("## Key Findings\n\n")
        f.write(f"- Total parameter groups: {len(param_groups)}\n")
        f.write(f"- Reinflation events: {len(reinflate_iters)}\n")

        if len(reinflate_iters) > 0:
            f.write(f"- Mean CV change at reinflation: ")
            cv_changes_at_reinflate = []
            for r_iter in reinflate_iters:
                before = evolution_df[evolution_df['Iteration'] == r_iter - 1]['CV'].mean()
                after = evolution_df[evolution_df['Iteration'] == r_iter]['CV'].mean()
                cv_changes_at_reinflate.append(after - before)
            f.write(f"{np.mean(cv_changes_at_reinflate):.2f}%\n")

        f.write("\n")

    print(f"Markdown summary saved to: {md_path}")

    print("\n" + "="*100)
    print("ANALYSIS COMPLETE!")
    print("="*100)

    return evolution_df


if __name__ == '__main__':
    # Get model name from command line or use default
    if len(sys.argv) > 1:
        model_name = sys.argv[1]
    else:
        model_name = DEFAULT_MODEL_NAME
        print(f"No model name provided, using MODEL_NAME from setup.py: {model_name}")

    analyze_reinflation_impact(model_name)
