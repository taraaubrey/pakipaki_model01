"""
s_recession_distributions.py - Plot recession distributions for prior vs posterior

Creates histogram plots for recession observations (pk4 and fliptime)
showing prior and posterior distributions for each stress period.

Usage:
    python scripts/s_recession_distributions.py run_name [options]

Examples:
    python scripts/s_recession_distributions.py local_run34 --post-iter 19
    python scripts/s_recession_distributions.py local_run34 --post-iter 19 --filter-file models/local_run34/filtered_realizations.csv
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

# Import default model name from setup
try:
    from setup import MODEL_NAME as DEFAULT_MODEL_NAME
except ImportError:
    DEFAULT_MODEL_NAME = None


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Plot recession distributions for prior vs posterior'
    )
    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default=DEFAULT_MODEL_NAME,
        help=f'Model name for file prefixes (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )

    parser.add_argument('--run-name', '-r', type=str, default=None,
                       help='Run name for directory path (default: same as model_name)')
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


def get_recession_obs(obs_data):
    """
    Get recession observation names organized by type and kper.

    Returns dict: {recession_type: {kper: obs_name}}
    """
    # Filter recession observations
    recession_obs = obs_data[obs_data['oname'] == 'recession'].copy()

    if len(recession_obs) == 0:
        print("  No recession observations found with oname='recession'")
        return {}

    recession_obs['kper'] = pd.to_numeric(recession_obs['kper'], errors='coerce')

    # Define recession types to look for in usecol field
    recession_types = {
        'pk4': ['pk4'],
        'fliptime': ['fliptime']
    }

    result = {}

    for rtype, usecol_patterns in recession_types.items():
        result[rtype] = {}

        for idx, row in recession_obs.iterrows():
            usecol = str(row.get('usecol', '')).lower()

            # Check if this observation matches the recession type
            for pattern in usecol_patterns:
                if pattern == usecol:
                    kper = int(row['kper']) if pd.notna(row['kper']) else 1
                    result[rtype][kper] = idx
                    break

    return result


def load_truth_recessions(model_ws):
    """
    Load truth recession values.

    Returns dict: {recession_type: {kper: value}}
    """
    # Truth file is in project root, not model directory
    truth_file = os.path.join('truth', 'output.sample_recession_rates.truth.csv')

    if not os.path.exists(truth_file):
        print(f"  Truth file not found: {truth_file}")
        return {}

    try:
        truth_df = pd.read_csv(truth_file)

        # Map column names to our recession type names
        col_mapping = {
            'pk4': 'pk4',
            'fliptime': 'fliptime'
        }

        # Convert to dict {recession_type: {kper: value}}
        truth_dict = {}
        for col, rtype in col_mapping.items():
            if col in truth_df.columns:
                truth_dict[rtype] = {}
                for _, row in truth_df.iterrows():
                    kper = int(row['kper'])
                    val = row[col]
                    if pd.notna(val):  # Include all non-NaN values
                        truth_dict[rtype][kper] = float(val)

        return truth_dict
    except Exception as e:
        print(f"  Error loading truth file: {e}")
        return {}


def create_recession_plots(run_name, model_name, post_iter=19, filter_file=None, suffix=None):
    """
    Create recession distribution plots comparing prior vs posterior.
    """
    # Set suffix for output files
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"RECESSION DISTRIBUTION ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    model_ws = os.path.join('models', run_name)
    master_dir = os.path.join(model_ws, 'pest', 'master_ies')
    output_dir = os.path.join(model_ws, 'pp_recession_analysis')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Get recession observations
    recession_obs = get_recession_obs(obs_data)

    print("\nRecession observations found:")
    all_obs_names = []
    for rtype, kper_dict in recession_obs.items():
        if kper_dict:
            print(f"  {rtype}: {kper_dict}")
            all_obs_names.extend(kper_dict.values())

    if len(all_obs_names) == 0:
        print("Error: No recession observations found")
        return

    # Load truth values
    print("\nLoading truth values...")
    truth_values = load_truth_recessions(model_ws)
    print(f"  Found {len(truth_values)} recession types with truth values")

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

    # Create plots for each recession type
    print("\nCreating plots...")

    recession_types_with_data = [rt for rt, kd in recession_obs.items() if len(kd) > 0]

    for rtype in recession_types_with_data:
        kper_dict = recession_obs[rtype]
        kpers = sorted(kper_dict.keys())

        if len(kpers) == 0:
            continue

        # Create figure with subplots for each kper (4 in a row)
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

            # Remove NaN values
            prior_vals = prior_vals[~np.isnan(prior_vals)]
            post_vals = post_vals[~np.isnan(post_vals)]

            if len(post_vals) == 0:
                ax.text(0.5, 0.5, 'No valid data', ha='center', va='center')
                ax.set_title(f'kper {kper}')
                continue

            # Separate bins for prior and posterior
            bins_prior = np.linspace(np.min(prior_vals), np.max(prior_vals), 25)
            bins_post = np.linspace(np.min(post_vals), np.max(post_vals), 25)

            # Plot histograms
            ax.hist(prior_vals, bins=bins_prior, alpha=0.5, color='#1f77b4',
                   label='Prior', edgecolor='black', linewidth=0.3)
            ax.hist(post_vals, bins=bins_post, alpha=0.5, color='#ff7f0e',
                   label='Posterior', edgecolor='black', linewidth=0.3)

            # Set x-axis limits based on posterior distribution
            post_min = np.min(post_vals)
            post_max = np.max(post_vals)
            x_margin = (post_max - post_min) * 0.1
            ax.set_xlim(post_min - x_margin, post_max + x_margin)

            # Add median lines
            ymin, ymax = ax.get_ylim()

            prior_median = np.median(prior_vals)
            post_median = np.median(post_vals)

            ax.axvline(prior_median, color='#1f77b4', linestyle='--', linewidth=2,
                      label=f'Prior median: {prior_median:.3f}')
            ax.axvline(post_median, color='#ff7f0e', linestyle='-', linewidth=2,
                      label=f'Post median: {post_median:.3f}')

            # Add truth value if available
            if rtype in truth_values and kper in truth_values[rtype]:
                truth_val = truth_values[rtype][kper]
                ax.axvline(truth_val, color='red', linestyle=':', linewidth=2,
                          label=f'Truth: {truth_val:.3f}')

            # Set axis labels based on type
            if rtype == 'pk4':
                ax.set_xlabel('Recession Rate (m/d)', fontsize=9)
            elif rtype == 'fliptime':
                ax.set_xlabel('Flip Time (days)', fontsize=9)
            else:
                ax.set_xlabel('Value', fontsize=9)

            ax.set_ylabel('Frequency', fontsize=9)
            ax.set_title(f'kper {kper}', fontsize=10, fontweight='bold')
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=7, loc='best')
            ax.tick_params(labelsize=8)

        fig.suptitle(f'{rtype} Recession', fontsize=12, fontweight='bold', y=1.02)
        plt.tight_layout()

        # Save figure
        output_path = os.path.join(output_dir, f'{run_name}_{rtype}_recession{file_suffix}.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {output_path}")

    print(f"\nAll plots saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_recession_plots(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        suffix=args.suffix
    )
