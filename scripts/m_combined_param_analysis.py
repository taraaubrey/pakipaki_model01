"""
Combined parameter analysis plot with evolution across iterations.

Combines elements from h_analyze_params.py and i_analyze_param_reinflation.py:
- Row 1: A - Lineplot of parameter change magnitude per iteration
- Row 2: B, C, D, E - Evolution of parameter mean for 4 main groups
- Row 3: F, G, H, I - Evolution of parameter CV for 4 main groups
- Row 4: J, K - Histograms of CV reduction (from prior and from last inflation)

Usage:
    python scripts/m_combined_param_analysis.py [model_name] [--run-name NAME] [--last-iter N]

    If model_name is not provided, uses MODEL_NAME from setup.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob
import pyemu
import argparse

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

# Import setup and parameter categorization
from setup import MODEL_NAME as DEFAULT_MODEL_NAME
from setup import REINFLATE_ITERS as DEFAULT_REINFLATE_LIST


def categorize_parameter(pname):
    """
    Categorize parameter by type based on pname prefix.

    Args:
        pname: Parameter name (e.g., 'npfklayer1-gr')

    Returns:
        tuple: (category, subcategory)
    """
    pname_lower = str(pname).lower()

    if 'npfk' in pname_lower:
        if '-gr' in pname_lower:
            return ('Hydraulic Conductivity', 'K-grid')
        elif '-pp' in pname_lower:
            return ('Hydraulic Conductivity', 'K-pilot')
        elif '-cn' in pname_lower:
            return ('Hydraulic Conductivity', 'K-constant')
        else:
            return ('Hydraulic Conductivity', 'K-other')

    elif 'stoss' in pname_lower or 'sto' in pname_lower:
        if '-gr' in pname_lower:
            return ('Storage', 'S-grid')
        elif '-pp' in pname_lower:
            return ('Storage', 'S-pilot')
        elif '-cn' in pname_lower:
            return ('Storage', 'S-constant')
        else:
            return ('Storage', 'S-other')

    elif 'rch' in pname_lower:
        if 'rchss' in pname_lower or 'ss' in pname_lower:
            return ('Recharge', 'RCH-steady')
        elif 'rchtr' in pname_lower or 'tr' in pname_lower:
            return ('Recharge', 'RCH-transient')
        else:
            return ('Recharge', 'RCH-other')

    elif 'ghb' in pname_lower:
        # Determine GHB type
        ghb_type = 'GHB-other'
        if 'ghbaw' in pname_lower:
            ghb_type = 'GHB-Awanui'
        elif 'ghbpw' in pname_lower:
            ghb_type = 'GHB-Poukawa'
        elif 'ghbspring' in pname_lower:
            ghb_type = 'GHB-Spring'
        elif 'ghbconf' in pname_lower:
            ghb_type = 'GHB-Confined'

        # Determine property
        if '-cond' in pname_lower:
            return ('GHB Conductance', ghb_type + '-cond')
        elif '-head' in pname_lower:
            return ('GHB Head', ghb_type + '-head')
        else:
            return ('GHB Other', ghb_type)

    else:
        return ('Other', 'Other')


def calculate_reinflation_iters(reinflate_list):
    """
    Calculate cumulative reinflation iteration numbers.

    Args:
        reinflate_list: List of iteration counts between reinflations

    Returns:
        list: Cumulative iteration numbers where reinflation occurred
    """
    r_iters = []
    start = 0
    for i in reinflate_list:
        r_iters.append(start + i + 1)
        start = start + i + 1
    return r_iters


def load_parameter_stats_all_iters(master_dir, model_name, pst):
    """
    Load parameter ensemble statistics for all iterations.

    Args:
        master_dir: Path to master_ies directory
        model_name: Model name
        pst: pyemu.Pst object

    Returns:
        DataFrame with iteration, group, category, mean, std, cv columns
    """
    # Find all parameter files
    par_files = sorted(glob.glob(os.path.join(master_dir, f'{model_name}.*.par.csv')))

    # Extract iteration numbers
    iterations = []
    for pf in par_files:
        basename = os.path.basename(pf)
        iter_num = int(basename.replace('.par.csv', '').split('.')[-1])
        iterations.append(iter_num)

    iterations = sorted(set(iterations))  # Remove duplicates
    print(f"Found parameter files for iterations: {iterations}")

    # Load parameter data
    param_data = pst.parameter_data.copy()

    # Create mapping from pargp to pname
    pargp_to_pname = {}
    for pargp in pst.par_groups:
        pars = param_data[param_data['pargp'] == pargp]
        if len(pars) > 0 and 'pname' in pars.columns:
            pargp_to_pname[pargp] = pars['pname'].iloc[0]
        else:
            pargp_to_pname[pargp] = pargp

    # Load statistics for all iterations
    all_stats = []

    for iter_num in iterations:
        par_file = os.path.join(master_dir, f'{model_name}.{iter_num}.par.csv')
        if not os.path.exists(par_file):
            continue

        # Load parameter ensemble
        par_df = pd.read_csv(par_file, index_col=0, low_memory=False)

        # Calculate statistics by parameter group
        for pargp in param_data['pargp'].unique():
            # Get all parameters in this group
            group_params = param_data[param_data['pargp'] == pargp].index.tolist()
            group_params = [p for p in group_params if p in par_df.columns]

            if len(group_params) == 0:
                continue

            # Get values for this group across all realizations
            group_values = par_df[group_params].values.flatten()

            # Calculate statistics
            mean_val = np.mean(group_values)
            std_val = np.std(group_values)
            cv_val = (std_val / mean_val * 100) if mean_val != 0 else np.nan

            # Get pname and category
            pname = pargp_to_pname.get(pargp, pargp)
            category, subcategory = categorize_parameter(pname)

            all_stats.append({
                'Iteration': iter_num,
                'Group': pargp,
                'Pname': pname,
                'Category': category,
                'Subcategory': subcategory,
                'Mean': mean_val,
                'Std': std_val,
                'CV': cv_val,
                'N_Params': len(group_params)
            })

    return pd.DataFrame(all_stats), iterations


def load_pcs_data(master_dir, model_name):
    """
    Load PCS (parameter change summary) data for all iterations.

    Args:
        master_dir: Path to master_ies directory
        model_name: Model name

    Returns:
        DataFrame with PCS data, or None if not available
    """
    pcs_files = sorted(glob.glob(os.path.join(master_dir, f'{model_name}.*.pcs.csv')))

    if len(pcs_files) == 0:
        return None

    pcs_data = []
    for pcs_file in pcs_files:
        basename = os.path.basename(pcs_file)
        # Skip reinflate files
        if '.reinflate.' in basename:
            continue
        try:
            iter_num = int(basename.replace('.pcs.csv', '').split('.')[-1])
        except ValueError:
            continue

        pcs = pd.read_csv(pcs_file)
        pcs['iteration'] = iter_num
        pcs_data.append(pcs)

    if len(pcs_data) == 0:
        return None

    return pd.concat(pcs_data, ignore_index=True)


def create_combined_plot(run_name, model_name, last_iter=None):
    """
    Create combined parameter analysis plot.

    Args:
        run_name: Name for directory path
        model_name: Name of the model run for file prefixes
        last_iter: Last iteration to use (default: auto-detect)
    """
    print("=" * 100)
    print(f"COMBINED PARAMETER ANALYSIS: {run_name}: {model_name}")
    print("=" * 100)

    # Set paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_param_analysis')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

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

    # Calculate reinflation iterations
    reinflate_iters = calculate_reinflation_iters(DEFAULT_REINFLATE_LIST)
    print(f"\nReinflation iterations: {reinflate_iters}")

    # Load parameter statistics for all iterations
    print("\nLoading parameter ensemble data...")
    stats_df, iterations = load_parameter_stats_all_iters(master_dir, model_name, pst)

    if len(stats_df) == 0:
        print("ERROR: No parameter statistics loaded")
        return

    # Determine last iteration
    if last_iter is None:
        last_iter = max(iterations)
    print(f"Using last iteration: {last_iter}")

    # Find last inflation iteration (the last one before or at last_iter)
    last_inflation_iter = 0
    for r_iter in reinflate_iters:
        if r_iter <= last_iter:
            last_inflation_iter = r_iter
        else:
            break

    print(f"Last inflation iteration: {last_inflation_iter}")

    # Load PCS data for mean change magnitudes
    pcs_df = load_pcs_data(master_dir, model_name)

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

    # Select 4 main categories for individual plots
    main_categories = ['Hydraulic Conductivity', 'Storage', 'Recharge', 'GHB Conductance']
    # Filter to categories that exist in the data
    available_categories = stats_df['Category'].unique()
    main_categories = [c for c in main_categories if c in available_categories]
    # Fill remaining slots with other categories
    for cat in available_categories:
        if cat not in main_categories and len(main_categories) < 4:
            main_categories.append(cat)

    print(f"Main categories for plots: {main_categories}")

    # Get subcategories for each main category
    subcategory_colors = {}
    for category in main_categories:
        cat_subcats = stats_df[stats_df['Category'] == category]['Subcategory'].unique()
        # Assign distinct colors to subcategories
        n_subcats = len(cat_subcats)
        base_color = category_colors.get(category, '#95a5a6')
        # Create color variations
        import matplotlib.colors as mcolors
        base_rgb = mcolors.to_rgb(base_color)
        for j, subcat in enumerate(cat_subcats):
            # Vary lightness
            factor = 0.5 + 0.5 * (j / max(1, n_subcats - 1)) if n_subcats > 1 else 0.75
            subcategory_colors[subcat] = tuple(min(1, c * factor + (1 - factor) * 0.3) for c in base_rgb)

    # Create figure with new layout (3 rows for A, B-E, F-I)
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.3)

    # Row 1 (A): Lineplot of parameter change magnitude per iteration
    ax_a = fig.add_subplot(gs[0, :])

    if pcs_df is not None and len(pcs_df) > 0:
        # Map pargp to category using stats_df
        group_to_cat = dict(zip(stats_df['Group'], stats_df['Category']))

        # Calculate mean change by category for each iteration
        for category in available_categories:
            cat_pcs = pcs_df[pcs_df['group'].map(group_to_cat) == category]
            if len(cat_pcs) > 0:
                change_by_iter = cat_pcs.groupby('iteration')['mean_change'].mean()
                color = category_colors.get(category, '#95a5a6')
                ax_a.plot(change_by_iter.index, change_by_iter.values, 'o-', linewidth=2,
                         markersize=5, label=category, color=color, alpha=0.8)

        # Mark reinflation iterations
        for r_iter in reinflate_iters:
            if r_iter <= last_iter:
                ax_a.axvline(x=r_iter, color='red', linestyle='--', linewidth=1, alpha=0.5)

        ax_a.set_xlabel('Iteration', fontsize=10)
        ax_a.set_ylabel('Mean Change Magnitude', fontsize=10)
        ax_a.set_title('A. Parameter change magnitude per iteration', fontsize=11, fontweight='bold', loc='left')
        ax_a.legend(fontsize=7, loc='best', ncol=3)
        ax_a.grid(True, alpha=0.3)
        ax_a.set_xlim(0, last_iter)
        ax_a.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    else:
        ax_a.text(0.5, 0.5, 'No PCS data available', ha='center', va='center',
                transform=ax_a.transAxes, fontsize=10, color='gray')

    # Row 2 (B, C, D, E): Evolution of parameter mean for 4 main groups (by subcategory)
    plot_labels = ['B', 'C', 'D', 'E']
    for i, category in enumerate(main_categories[:4]):
        ax = fig.add_subplot(gs[1, i])
        cat_data = stats_df[stats_df['Category'] == category]

        # Plot mean evolution for each subcategory
        for subcat in cat_data['Subcategory'].unique():
            subcat_data = cat_data[cat_data['Subcategory'] == subcat]
            mean_by_iter = subcat_data.groupby('Iteration')['Mean'].mean()
            color = subcategory_colors.get(subcat, '#95a5a6')
            # Create shorter label
            label = subcat.replace(category + '-', '').replace('K-', '').replace('S-', '').replace('RCH-', '').replace('GHB-', '')
            ax.plot(mean_by_iter.index, mean_by_iter.values, 'o-', linewidth=1.5,
                   markersize=3, color=color, alpha=0.8, label=label)

        # Mark reinflation iterations
        for r_iter in reinflate_iters:
            if r_iter <= last_iter:
                ax.axvline(x=r_iter, color='red', linestyle='--', linewidth=1, alpha=0.3)

        ax.set_xlabel('Iteration', fontsize=8)
        if i == 0:
            ax.set_ylabel('Mean', fontsize=8)
        ax.set_title(f'{plot_labels[i]}. {category}', fontsize=9, fontweight='bold', loc='left')
        ax.legend(fontsize=5, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, last_iter)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.tick_params(labelsize=7)
        if category in ['Hydraulic Conductivity', 'GHB Conductance']:
            ax.set_yscale('log')

    # Row 3 (F, G, H, I): Evolution of parameter CV for 4 main groups (by subcategory)
    plot_labels = ['F', 'G', 'H', 'I']
    for i, category in enumerate(main_categories[:4]):
        ax = fig.add_subplot(gs[2, i])
        cat_data = stats_df[stats_df['Category'] == category]

        # Plot CV evolution for each subcategory
        for subcat in cat_data['Subcategory'].unique():
            subcat_data = cat_data[cat_data['Subcategory'] == subcat]
            cv_by_iter = subcat_data.groupby('Iteration')['CV'].mean()
            color = subcategory_colors.get(subcat, '#95a5a6')
            # Create shorter label
            label = subcat.replace(category + '-', '').replace('K-', '').replace('S-', '').replace('RCH-', '').replace('GHB-', '')
            ax.plot(cv_by_iter.index, cv_by_iter.values, 'o-', linewidth=1.5,
                   markersize=3, color=color, alpha=0.8, label=label)

        # Mark reinflation iterations
        for r_iter in reinflate_iters:
            if r_iter <= last_iter:
                ax.axvline(x=r_iter, color='red', linestyle='--', linewidth=1, alpha=0.3)

        ax.set_xlabel('Iteration', fontsize=8)
        if i == 0:
            ax.set_ylabel('CV (%)', fontsize=8)
        ax.set_title(f'{plot_labels[i]}. {category} CV', fontsize=9, fontweight='bold', loc='left')
        ax.legend(fontsize=5, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, last_iter)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
        ax.tick_params(labelsize=7)

    # Save first figure (A-I)
    fig_path = os.path.join(output_dir, f'{model_name}_combined_param_analysis.png')
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"\nFigure 1 saved to: {fig_path}")
    plt.close()

    # Create second figure for CV reduction bar charts (A and B)
    fig2 = plt.figure(figsize=(14, 8))
    gs2 = fig2.add_gridspec(1, 2, wspace=0.4)

    # Get CV at iteration 0 and last iteration for all groups
    prior_data = stats_df[stats_df['Iteration'] == 0].set_index('Group')
    final_data = stats_df[stats_df['Iteration'] == last_iter].set_index('Group')

    # Calculate CV reduction from prior for each subcategory
    cv_reduction_data = []
    for group in prior_data.index:
        if group in final_data.index:
            prior_cv = prior_data.loc[group, 'CV']
            final_cv = final_data.loc[group, 'CV']
            category = prior_data.loc[group, 'Category']
            # Skip GHB Head category
            if category == 'GHB Head':
                continue
            if prior_cv > 0:
                reduction = ((prior_cv - final_cv) / prior_cv * 100)
                subcat = prior_data.loc[group, 'Subcategory']

                # Create clearer label with pp/grid/constant
                if '-grid' in subcat.lower() or subcat.endswith('-gr'):
                    param_type = 'grid'
                elif '-pilot' in subcat.lower() or subcat.endswith('-pp'):
                    param_type = 'pp'
                elif '-constant' in subcat.lower() or subcat.endswith('-cn'):
                    param_type = 'const'
                else:
                    param_type = ''

                # Create display name
                base_name = subcat.replace('K-', '').replace('S-', '').replace('RCH-', '').replace('GHB-', '')
                base_name = base_name.replace('-grid', '').replace('-pilot', '').replace('-constant', '')
                base_name = base_name.replace('-gr', '').replace('-pp', '').replace('-cn', '')
                if param_type:
                    display_name = f"{base_name} ({param_type})"
                else:
                    display_name = base_name

                cv_reduction_data.append({
                    'Subcategory': subcat,
                    'DisplayName': display_name,
                    'Category': category,
                    'Reduction': reduction
                })

    if len(cv_reduction_data) > 0:
        cv_df = pd.DataFrame(cv_reduction_data)

        # Split into positive (reductions) and negative (increases)
        positive_df = cv_df[cv_df['Reduction'] >= 0].sort_values('Reduction', ascending=False)
        negative_df = cv_df[cv_df['Reduction'] < 0].sort_values('Reduction', ascending=True)

        # A: Positive reductions (CV decreased) - vertical bars
        ax_a = fig2.add_subplot(gs2[0, 0])

        if len(positive_df) > 0:
            colors = [subcategory_colors.get(subcat, category_colors.get(cat, '#95a5a6'))
                      for subcat, cat in zip(positive_df['Subcategory'], positive_df['Category'])]

            x_pos = range(len(positive_df))
            ax_a.bar(x_pos, positive_df['Reduction'].values, color=colors, alpha=0.8, edgecolor='black')

            ax_a.set_xticks(x_pos)
            ax_a.set_xticklabels(positive_df['DisplayName'].values, rotation=45, ha='right', fontsize=8)

            ax_a.set_ylabel('CV Reduction (%)', fontsize=10)
            ax_a.set_title(f'A. CV reductions (iter 0 → {last_iter})', fontsize=11, fontweight='bold', loc='left')
            ax_a.axhline(y=0, color='gray', linestyle=':', linewidth=1)
            ax_a.grid(True, alpha=0.3, axis='y')
            ax_a.set_axisbelow(True)
        else:
            ax_a.text(0.5, 0.5, 'No CV reductions', ha='center', va='center',
                     transform=ax_a.transAxes, fontsize=10, color='gray')

        # B: Negative reductions (CV increased) - vertical bars
        ax_b = fig2.add_subplot(gs2[0, 1])

        if len(negative_df) > 0:
            colors = [subcategory_colors.get(subcat, category_colors.get(cat, '#95a5a6'))
                      for subcat, cat in zip(negative_df['Subcategory'], negative_df['Category'])]

            x_pos = range(len(negative_df))
            ax_b.bar(x_pos, negative_df['Reduction'].values, color=colors, alpha=0.8, edgecolor='black')

            ax_b.set_xticks(x_pos)
            ax_b.set_xticklabels(negative_df['DisplayName'].values, rotation=45, ha='right', fontsize=8)

            ax_b.set_ylabel('CV Reduction (%)', fontsize=10)
            ax_b.set_title(f'B. CV increases (iter 0 → {last_iter})', fontsize=11, fontweight='bold', loc='left')
            ax_b.axhline(y=0, color='gray', linestyle=':', linewidth=1)
            ax_b.grid(True, alpha=0.3, axis='y')
            ax_b.set_axisbelow(True)
        else:
            ax_b.text(0.5, 0.5, 'No CV increases', ha='center', va='center',
                     transform=ax_b.transAxes, fontsize=10, color='gray')

    # Save second figure
    fig2_path = os.path.join(output_dir, f'{model_name}_cv_reduction.png')
    plt.savefig(fig2_path, dpi=200, bbox_inches='tight')
    print(f"Figure 2 saved to: {fig2_path}")
    plt.close()

    # Create individual histogram figures for each parameter group
    # Compare final iteration vs reinflated prior (last inflation iteration)
    print("\nCreating parameter distribution histograms...")

    # Load parameter ensembles for final and reinflated iterations
    final_par_file = os.path.join(master_dir, f'{model_name}.{last_iter}.par.csv')
    prior_par_file = os.path.join(master_dir, f'{model_name}.{last_inflation_iter}.par.csv')

    if os.path.exists(final_par_file) and os.path.exists(prior_par_file):
        final_par = pd.read_csv(final_par_file, index_col=0, low_memory=False)
        prior_par = pd.read_csv(prior_par_file, index_col=0, low_memory=False)

        # Get parameter data
        param_data = pst.parameter_data.copy()

        # Create histogram subdirectory
        hist_dir = os.path.join(output_dir, 'histograms')
        if not os.path.exists(hist_dir):
            os.makedirs(hist_dir)
            print(f"Created histogram directory: {hist_dir}")

        # Get unique subcategories
        subcategories = stats_df[['Subcategory', 'Category']].drop_duplicates()

        for _, row in subcategories.iterrows():
            subcat = row['Subcategory']
            category = row['Category']

            # Skip GHB Head
            if category == 'GHB Head':
                continue

            # Get parameter groups for this subcategory
            subcat_groups = stats_df[stats_df['Subcategory'] == subcat]['Group'].unique()

            # Collect all parameter values for this subcategory
            final_values = []
            prior_values = []

            for group in subcat_groups:
                # Get parameters in this group
                group_params = param_data[param_data['pargp'] == group].index.tolist()
                group_params = [p for p in group_params if p in final_par.columns and p in prior_par.columns]

                if len(group_params) > 0:
                    final_values.extend(final_par[group_params].values.flatten())
                    prior_values.extend(prior_par[group_params].values.flatten())

            if len(final_values) > 0 and len(prior_values) > 0:
                # Create figure (1/3 page size ~ 4x3 inches)
                fig3, ax = plt.subplots(figsize=(4, 3))

                # Determine if log scale needed
                all_vals = np.concatenate([final_values, prior_values])
                use_log = np.min(all_vals[all_vals > 0]) > 0 and np.max(all_vals) / np.min(all_vals[all_vals > 0]) > 100

                if use_log:
                    # Log-transform for histogram
                    final_log = np.log10(np.array(final_values)[np.array(final_values) > 0])
                    prior_log = np.log10(np.array(prior_values)[np.array(prior_values) > 0])

                    # Create bins
                    all_log = np.concatenate([final_log, prior_log])
                    bins = np.linspace(np.min(all_log), np.max(all_log), 30)

                    # Plot histograms
                    ax.hist(prior_log, bins=bins, alpha=0.5, color='blue',
                           label=f'Iter {last_inflation_iter} (reinflated)', edgecolor='black', linewidth=0.5)
                    ax.hist(final_log, bins=bins, alpha=0.5, color='red',
                           label=f'Iter {last_iter} (final)', edgecolor='black', linewidth=0.5)

                    ax.set_xlabel('log₁₀(value)', fontsize=9)
                else:
                    # Use linear scale
                    bins = 30

                    ax.hist(prior_values, bins=bins, alpha=0.5, color='blue',
                           label=f'Iter {last_inflation_iter} (reinflated)', edgecolor='black', linewidth=0.5)
                    ax.hist(final_values, bins=bins, alpha=0.5, color='red',
                           label=f'Iter {last_iter} (final)', edgecolor='black', linewidth=0.5)

                    ax.set_xlabel('Value', fontsize=9)

                ax.set_ylabel('Frequency', fontsize=9)
                ax.set_title(f'{subcat}', fontsize=10, fontweight='bold')
                ax.legend(fontsize=7, loc='best')
                ax.grid(True, alpha=0.3)
                ax.tick_params(labelsize=8)

                # Save figure
                safe_name = subcat.replace('/', '_').replace('\\', '_').replace(' ', '_')
                hist_path = os.path.join(hist_dir, f'{safe_name}_hist.png')
                plt.savefig(hist_path, dpi=150, bbox_inches='tight')
                plt.close()

        print(f"Histogram figures saved to: {hist_dir}")
    else:
        print(f"WARNING: Could not find parameter files for iterations {last_iter} and {last_inflation_iter}")

    # Save data to CSV
    csv_path = os.path.join(output_dir, f'{model_name}_param_stats_all_iters.csv')
    stats_df.to_csv(csv_path, index=False)
    print(f"Statistics saved to: {csv_path}")

    print("\n" + "=" * 100)
    print("ANALYSIS COMPLETE!")
    print("=" * 100)

    return stats_df


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Combined parameter analysis plot',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default=DEFAULT_MODEL_NAME,
        help=f'Model name for file prefixes (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )

    parser.add_argument(
        '--run-name', '-r',
        type=str,
        default=None,
        help='Run name for directory path (default: same as model_name)'
    )

    parser.add_argument(
        '--last-iter', '-i',
        type=int,
        default=None,
        help='Last iteration to use (default: auto-detect)'
    )

    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_combined_plot(run_name, model_name, args.last_iter)
