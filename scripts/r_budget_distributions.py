"""
r_budget_distributions.py - Plot budget distributions for prior vs posterior

Creates histogram plots for each budget type (Awanui, Poukawa, Spring, Confined,
Recharge, Storage) showing prior and posterior distributions for each stress period.

Usage:
    python scripts/r_budget_distributions.py run_name [options]

Examples:
    python scripts/r_budget_distributions.py local_run34 --post-iter 19
    python scripts/r_budget_distributions.py local_run34 --post-iter 19 --filter-file models/local_run34/filtered_realizations.csv
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


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Plot budget distributions for prior vs posterior'
    )
    parser.add_argument('run_name', type=str, help='Name of the run directory')
    
    parser.add_argument('--model-name', type=str, default=None,
                       help='Model name (default: same as run_name)')
    parser.add_argument('--post-iter', type=int, default=19,
                       help='Posterior iteration (default: 19)')
    parser.add_argument('--filter-file', type=str, default=None,
                       help='Path to file with filtered realization indices')
    parser.add_argument('--suffix', type=str, default=None,
                       help='Suffix for output files (e.g., "_filtered")')
    return parser.parse_args()


def load_filtered_realizations(filter_file):
    """Load filtered realization indices from file."""
    if filter_file is None or not os.path.exists(filter_file):
        return None
    df = pd.read_csv(filter_file, comment='#')
    return df['realization'].tolist()


def get_budget_obs(obs_data):
    """
    Get budget observation names organized by type and kper.

    Returns dict: {budget_type: {kper: obs_name}}
    """
    # Filter budget observations
    budget_obs = obs_data[obs_data['oname'] == 'budget'].copy()

    if len(budget_obs) == 0:
        print("  No budget observations found with oname='budget'")
        return {}

    budget_obs['kper'] = pd.to_numeric(budget_obs['kper'], errors='coerce')

    # Define budget types to look for in usecol field
    budget_types = {
        'Awanui': ['awanui'],
        'Poukawa': ['poukawa'],
        'Spring': ['spring'],
        'Confined': ['confined'],
        'Recharge': ['recharge'],
        'Storage': ['stoss']
    }

    result = {}

    for btype, usecol_patterns in budget_types.items():
        result[btype] = {}

        for idx, row in budget_obs.iterrows():
            usecol = str(row.get('usecol', '')).lower()

            # Check if this observation matches the budget type
            for pattern in usecol_patterns:
                if pattern == usecol:
                    kper = int(row['kper']) if pd.notna(row['kper']) else 1
                    result[btype][kper] = idx
                    break

    return result


def load_truth_budgets(model_ws):
    """
    Load truth budget values.

    Returns dict: {budget_type: {kper: value}}
    """
    # Truth file is in project root, not model directory
    truth_file = os.path.join('truth', 'output.budget.truth.csv')

    if not os.path.exists(truth_file):
        print(f"  Truth file not found: {truth_file}")
        return {}

    try:
        truth_df = pd.read_csv(truth_file)

        # Map column names to our budget type names
        col_mapping = {
            'awanui': 'Awanui',
            'poukawa': 'Poukawa',
            'spring': 'Spring',
            'confined': 'Confined',
            'recharge': 'Recharge',
            'stoss': 'Storage'
        }

        # Convert to dict {budget_type: {kper: value}}
        truth_dict = {}
        for col, btype in col_mapping.items():
            if col in truth_df.columns:
                truth_dict[btype] = {}
                for _, row in truth_df.iterrows():
                    kper = int(row['kper'])
                    val = row[col]
                    if val != 0:  # Only include non-zero values
                        truth_dict[btype][kper] = float(val)

        return truth_dict
    except Exception as e:
        print(f"  Error loading truth file: {e}")
        return {}


def create_budget_plots(run_name, model_name, post_iter=19, filter_file=None, suffix=None):
    """
    Create budget distribution plots comparing prior vs posterior.
    """
    # Set suffix for output files
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"BUDGET DISTRIBUTION ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    model_ws = os.path.join('models', run_name)
    master_dir = os.path.join(model_ws, 'pest', 'master_ies')
    output_dir = os.path.join(model_ws, 'pp_budget_analysis')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Get budget observations
    budget_obs = get_budget_obs(obs_data)

    print("\nBudget observations found:")
    all_obs_names = []
    for btype, kper_dict in budget_obs.items():
        if kper_dict:
            print(f"  {btype}: {kper_dict}")
            all_obs_names.extend(kper_dict.values())

    if len(all_obs_names) == 0:
        print("Error: No budget observations found")
        return

    # Load truth values
    print("\nLoading truth values...")
    truth_values = load_truth_budgets(model_ws)
    print(f"  Found {len(truth_values)} non-zero truth values")

    # Load filtered realizations
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations")

    # Load observation ensembles
    print(f"\nLoading observation ensembles...")
    prior_file = os.path.join(master_dir, f'{model_name}.0.obs.csv')
    post_file = os.path.join(master_dir, f'{model_name}.{post_iter}.obs.csv')

    # Read only needed columns
    print(f"  Reading prior (iter 0)...")
    prior_header = pd.read_csv(prior_file, nrows=0)
    cols_to_read = [prior_header.columns[0]]
    cols_to_read.extend([col for col in all_obs_names if col in prior_header.columns])
    prior_oe = pd.read_csv(prior_file, usecols=cols_to_read, index_col=0)

    print(f"  Reading posterior (iter {post_iter})...")
    post_oe = pd.read_csv(post_file, usecols=cols_to_read, index_col=0)

    # Filter realizations
    if filtered_ids:
        prior_oe = prior_oe.loc[[r for r in filtered_ids if r in prior_oe.index]]
        post_oe = post_oe.loc[[r for r in filtered_ids if r in post_oe.index]]
        print(f"  Filtered to {len(prior_oe)} prior, {len(post_oe)} posterior realizations")

    # Save extracted data for quick reloading
    print("\nSaving extracted data...")
    data_dir = os.path.join(output_dir, 'data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Build dataframe with all budget data
    budget_data_rows = []
    for btype in budget_types_with_data:
        kper_dict = budget_obs[btype]
        for kper, obs_name in kper_dict.items():
            if obs_name in prior_oe.columns:
                for real_id in prior_oe.index:
                    budget_data_rows.append({
                        'budget_type': btype,
                        'kper': kper,
                        'realization': real_id,
                        'iteration': 'prior',
                        'value': prior_oe.loc[real_id, obs_name]
                    })
                for real_id in post_oe.index:
                    budget_data_rows.append({
                        'budget_type': btype,
                        'kper': kper,
                        'realization': real_id,
                        'iteration': 'posterior',
                        'value': post_oe.loc[real_id, obs_name]
                    })

    budget_data_df = pd.DataFrame(budget_data_rows)
    budget_data_file = os.path.join(data_dir, f'{run_name}_budget_data{file_suffix}.csv')
    budget_data_df.to_csv(budget_data_file, index=False)
    print(f"  Saved: {budget_data_file}")

    # Save truth values
    if truth_values:
        truth_rows = []
        for btype, kper_vals in truth_values.items():
            for kper, val in kper_vals.items():
                truth_rows.append({
                    'budget_type': btype,
                    'kper': kper,
                    'truth_value': val
                })
        truth_df = pd.DataFrame(truth_rows)
        truth_data_file = os.path.join(data_dir, f'{run_name}_budget_truth{file_suffix}.csv')
        truth_df.to_csv(truth_data_file, index=False)
        print(f"  Saved: {truth_data_file}")

    # Create plots for each budget type
    print("\nCreating plots...")

    budget_types_with_data = [bt for bt, kd in budget_obs.items() if len(kd) > 0]

    for btype in budget_types_with_data:
        kper_dict = budget_obs[btype]
        kpers = sorted(kper_dict.keys())

        if len(kpers) == 0:
            continue

        # Create figure with subplots for each kper
        n_kpers = len(kpers)
        fig, axes = plt.subplots(1, n_kpers, figsize=(4*n_kpers, 4))

        if n_kpers == 1:
            axes = [axes]

        for ax, kper in zip(axes, kpers):
            obs_name = kper_dict[kper]

            if obs_name not in prior_oe.columns:
                ax.text(0.5, 0.5, 'Data not available', ha='center', va='center')
                ax.set_title(f'kper {kper}')
                continue

            # Get data
            prior_vals = prior_oe[obs_name].values
            post_vals = post_oe[obs_name].values

            # Determine bins
            all_vals = np.concatenate([prior_vals, post_vals])
            bins = np.linspace(np.min(all_vals), np.max(all_vals), 25)

            # Plot histograms
            ax.hist(prior_vals, bins=bins, alpha=0.5, color='#1f77b4',
                   label='Prior', edgecolor='black', linewidth=0.3)
            ax.hist(post_vals, bins=bins, alpha=0.5, color='#ff7f0e',
                   label='Posterior', edgecolor='black', linewidth=0.3)

            # Add mean lines
            ymin, ymax = ax.get_ylim()

            prior_mean = np.mean(prior_vals)
            post_mean = np.mean(post_vals)

            ax.axvline(prior_mean, color='#1f77b4', linestyle='--', linewidth=2,
                      label=f'Prior mean: {prior_mean:.1f}')
            ax.axvline(post_mean, color='#ff7f0e', linestyle='-', linewidth=2,
                      label=f'Post mean: {post_mean:.1f}')

            # Add truth value if available
            if btype in truth_values and kper in truth_values[btype]:
                truth_val = truth_values[btype][kper]
                ax.axvline(truth_val, color='red', linestyle=':', linewidth=2,
                          label=f'Truth: {truth_val:.1f}')

            ax.set_xlabel('Budget (m³/d)', fontsize=9)
            ax.set_ylabel('Frequency', fontsize=9)
            ax.set_title(f'kper {kper}', fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7, loc='best')
            ax.tick_params(labelsize=8)

        fig.suptitle(f'{btype} Budget', fontsize=12, fontweight='bold', y=1.02)
        plt.tight_layout()

        # Save figure
        output_path = os.path.join(output_dir, f'{run_name}_{btype.lower()}_budget{file_suffix}.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {output_path}")

    print(f"\nAll plots saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    run_name = args.run_name
    model_name = args.model_name if args.model_name else run_name

    create_budget_plots(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        suffix=args.suffix
    )
