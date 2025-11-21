"""
n_param_histograms.py - Generate parameter distribution histograms

Creates individual histogram figures for each parameter group/subgroup comparing:
- Final iteration distribution
- Reinflated prior distribution (last inflation iteration)

Usage:
    python scripts/n_param_histograms.py run_name [--last-iter N]

Example:
    python scripts/n_param_histograms.py local_run34 --last-iter 19
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemu

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from setup import REINFLATE_ITERS
from h_analyze_params import categorize_parameter


def calculate_reinflation_iters(reinflate_list):
    """Calculate iteration indices where reinflation occurred."""
    reinflation_iters = []
    cumulative = 0
    for n_iters in reinflate_list:
        cumulative += n_iters
        reinflation_iters.append(cumulative)
    return reinflation_iters


def create_param_histograms(run_name, model_name, last_iter=None, prior_iter=None,
                            filter_file=None, suffix=None):
    """
    Create individual histogram figures for each parameter subcategory.

    Parameters:
    -----------
    run_name : str
        Name of the run directory
    model_name : str
        Name of the model
    last_iter : int, optional
        Last iteration index (default: auto-detect)
    prior_iter : int, optional
        Prior/reinflated iteration to compare against (default: auto-calculate)
    filter_file : str, optional
        Path to file with filtered realization indices
    suffix : str, optional
        Suffix for output files (e.g., "_filtered")
    """

    # Set suffix for output files
    file_suffix = suffix if suffix else ''

    # Set up paths
    master_dir = f'models/{run_name}/pest/master_ies'
    output_dir = f'models/{run_name}/pp_param_analysis'
    hist_dir = os.path.join(output_dir, 'histograms')

    # Create output directories
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(hist_dir):
        os.makedirs(hist_dir)

    # Load PST file
    pst_path = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_path)
    param_data = pst.parameter_data.copy()

    # Calculate reinflation iterations (used if prior_iter not specified)
    reinflation_iters = calculate_reinflation_iters(REINFLATE_ITERS)

    # Use provided prior_iter or calculate from reinflation iterations
    if prior_iter is not None:
        last_inflation_iter = prior_iter
    else:
        last_inflation_iter = reinflation_iters[-1] if reinflation_iters else 0

    # Auto-detect last iteration if not specified
    if last_iter is None:
        iter_files = [f for f in os.listdir(master_dir) if f.endswith('.par.csv')]
        iters = []
        for f in iter_files:
            try:
                iter_num = int(f.split('.')[1])
                iters.append(iter_num)
            except (ValueError, IndexError):
                pass
        last_iter = max(iters) if iters else 19

    print(f"Run: {run_name}")
    print(f"Model: {model_name}")
    print(f"Final iteration: {last_iter}")
    print(f"Last inflation iteration: {last_inflation_iter}")
    print(f"Output directory: {hist_dir}")

    # Load parameter ensembles
    final_par_file = os.path.join(master_dir, f'{model_name}.{last_iter}.par.csv')
    prior_par_file = os.path.join(master_dir, f'{model_name}.{last_inflation_iter}.par.csv')

    if not os.path.exists(final_par_file):
        print(f"Error: Final parameter file not found: {final_par_file}")
        return
    if not os.path.exists(prior_par_file):
        print(f"Error: Prior parameter file not found: {prior_par_file}")
        return

    print("\nLoading parameter ensembles...")
    final_par = pd.read_csv(final_par_file, index_col=0, low_memory=False)
    prior_par = pd.read_csv(prior_par_file, index_col=0, low_memory=False)

    # Load and apply filtered realizations if specified
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        final_par = final_par.loc[[r for r in filtered_ids if r in final_par.index]]
        prior_par = prior_par.loc[[r for r in filtered_ids if r in prior_par.index]]
        print(f"Filtered to {len(final_par)} final, {len(prior_par)} prior realizations")

    # Create mapping from pargp to pname (using first parameter name in each group)
    pargp_to_pname = {}
    for pargp in param_data['pargp'].unique():
        pars = param_data[param_data['pargp'] == pargp]
        if len(pars) > 0 and 'pname' in pars.columns:
            pargp_to_pname[pargp] = pars['pname'].iloc[0]
        else:
            pargp_to_pname[pargp] = pargp

    # Build subcategory mapping
    subcat_groups = {}
    for pargp in param_data['pargp'].unique():
        pname = pargp_to_pname.get(pargp, pargp)
        category, subcat = categorize_parameter(pname)

        # Skip GHB Head
        if category == 'GHB Head':
            continue

        if subcat not in subcat_groups:
            subcat_groups[subcat] = {'category': category, 'groups': []}
        subcat_groups[subcat]['groups'].append(pargp)

    print(f"\nFound {len(subcat_groups)} parameter subcategories")
    print("\nCreating histograms...")

    # Define colors for categories
    category_colors = {
        'Hydraulic Conductivity': {'prior': "#8d8d8d", 'final': '#ff7f0e'},  # blue/orange
        'Storage': {'prior': '#8d8d8d', 'final': '#d62728'},  # green/red
        'Recharge': {'prior': '#8d8d8d', 'final': '#8c564b'},  # purple/brown
        'GHB Conductance': {'prior': '#8d8d8d', 'final': '#bcbd22'},  # cyan/olive
        'Other': {'prior': '#8d8d8d', 'final': '#e377c2'}  # gray/pink
    }

    created_count = 0

    for subcat, info in sorted(subcat_groups.items()):
        category = info['category']
        group_list = info['groups']

        # Collect all parameter values for this subcategory
        final_values = []
        prior_values = []

        for group in group_list:
            group_params = param_data[param_data['pargp'] == group].index.tolist()
            group_params = [p for p in group_params if p in final_par.columns and p in prior_par.columns]

            if len(group_params) > 0:
                final_values.extend(final_par[group_params].values.flatten())
                prior_values.extend(prior_par[group_params].values.flatten())

        if len(final_values) == 0 or len(prior_values) == 0:
            print(f"  Skipping {subcat}: no data")
            continue

        # Create figure (1/3 page size ~ 4x3 inches)
        fig, ax = plt.subplots(figsize=(4, 3))

        # Get colors for this category
        colors = category_colors.get(category, category_colors['Other'])

        # Determine if log scale needed
        all_vals = np.concatenate([final_values, prior_values])
        positive_vals = all_vals[all_vals > 0]

        if len(positive_vals) > 0:
            val_range = np.max(positive_vals) / np.min(positive_vals)
            use_log = val_range > 100
        else:
            use_log = False

        if use_log:
            # Log-transform for histogram
            final_positive = np.array(final_values)[np.array(final_values) > 0]
            prior_positive = np.array(prior_values)[np.array(prior_values) > 0]

            if len(final_positive) == 0 or len(prior_positive) == 0:
                print(f"  Skipping {subcat}: no positive values for log scale")
                plt.close()
                continue

            final_log = np.log10(final_positive)
            prior_log = np.log10(prior_positive)

            all_log = np.concatenate([final_log, prior_log])
            bins = np.linspace(np.min(all_log), np.max(all_log), 30)

            ax.hist(prior_log, bins=bins, alpha=0.5, color=colors['prior'],
                   label=f'Iter {last_inflation_iter} (reinflated)', edgecolor='black', linewidth=0.5)
            ax.hist(final_log, bins=bins, alpha=0.5, color=colors['final'],
                   label=f'Iter {last_iter} (final)', edgecolor='black', linewidth=0.5)

            ax.set_xlabel('log₁₀(value)', fontsize=9)
        else:
            bins = 30
            ax.hist(prior_values, bins=bins, alpha=0.5, color=colors['prior'],
                   label=f'Iter {last_inflation_iter} (reinflated)', edgecolor='black', linewidth=0.5)
            ax.hist(final_values, bins=bins, alpha=0.5, color=colors['final'],
                   label=f'Iter {last_iter} (final)', edgecolor='black', linewidth=0.5)
            ax.set_xlabel('Value', fontsize=9)

        ax.set_ylabel('Frequency', fontsize=9)
        ax.set_title(f'{subcat}', fontsize=10, fontweight='bold')
        ax.legend(fontsize=7, loc='best')
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)

        plt.tight_layout()

        # Save figure
        safe_name = subcat.replace('/', '_').replace('\\', '_').replace(' ', '_')
        hist_path = os.path.join(hist_dir, f'{run_name}_{safe_name}_hist{file_suffix}.png')
        plt.savefig(hist_path, dpi=150, bbox_inches='tight')
        plt.close()

        created_count += 1
        print(f"  Created: {safe_name}_hist{file_suffix}.png")

    print(f"\nCompleted: {created_count} histograms saved to {hist_dir}")


def load_filtered_realizations(filter_file):
    """Load filtered realization indices from file."""
    if filter_file is None or not os.path.exists(filter_file):
        return None
    df = pd.read_csv(filter_file, comment='#')
    return df['realization'].tolist()


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate parameter distribution histograms comparing final vs reinflated prior'
    )
    parser.add_argument('run_name', type=str, help='Name of the run directory')
    parser.add_argument('--model-name', type=str, default=None,
                       help='Model name (default: same as run_name)')
    parser.add_argument('--last-iter', type=int, default=None,
                       help='Last iteration index (default: auto-detect)')
    parser.add_argument('--prior-iter', type=int, default=None,
                       help='Prior/reinflated iteration to compare against (default: auto-calculate from REINFLATE_ITERS)')
    parser.add_argument('--filter-file', type=str, default=None,
                       help='Path to file with filtered realization indices')
    parser.add_argument('--suffix', type=str, default=None,
                       help='Suffix for output files (e.g., "_filtered")')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    run_name = args.run_name
    model_name = args.model_name if args.model_name else run_name

    create_param_histograms(
        run_name,
        model_name,
        last_iter=args.last_iter,
        prior_iter=args.prior_iter,
        filter_file=args.filter_file,
        suffix=args.suffix
    )
