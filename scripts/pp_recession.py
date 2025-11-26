"""
pp_recession.py - Compare prior vs posterior recession rate distributions

Compares prior (iter 0) and posterior distributions for recession rate predictions
showing present (kper 2) vs past (kper 4).

Usage:
    python scripts/pp_recession.py run_name [options]

Examples:
    python scripts/pp_recession.py local_run34 --post-iter 19
    python scripts/pp_recession.py run36 --run-name run36_62859_noGHBheads --post-iter 39
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
        description='Compare prior vs posterior recession rate distributions'
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
    parser.add_argument('--pct-low', type=int, default=10,
                       help='Lower percentile for filtering (default: 10)')
    parser.add_argument('--pct-high', type=int, default=90,
                       help='Upper percentile for filtering (default: 90)')
    return parser.parse_args()


def load_filtered_realizations(filter_file):
    """Load filtered realization indices from file."""
    if filter_file is None or not os.path.exists(filter_file):
        return None

    df = pd.read_csv(filter_file, comment='#')
    return df['realization'].tolist()


def get_recession_obs_by_kper(obs_data):
    """
    Get recession observation names organized by kper and usecol (site).

    Returns dict: {kper: {usecol: obs_name}}
    """
    recession_obs = obs_data[obs_data['oname'] == 'recession'].copy()
    recession_obs['kper'] = pd.to_numeric(recession_obs['kper'], errors='coerce')

    obs_dict = {}
    for idx, row in recession_obs.iterrows():
        kper = int(row['kper']) if pd.notna(row['kper']) else 1
        usecol = str(row.get('usecol', ''))

        if kper not in obs_dict:
            obs_dict[kper] = {}
        obs_dict[kper][usecol] = idx

    return obs_dict, recession_obs


def load_truth_data(run_name):
    """
    Load truth recession rate data.

    Returns dict: {kper: {usecol: value}}
    """
    truth_dir = os.path.join('models', run_name, 'truth')

    # Look for recession truth file
    truth_file = None
    if os.path.exists(truth_dir):
        for f in os.listdir(truth_dir):
            if 'recession' in f.lower() and f.endswith('.csv'):
                truth_file = os.path.join(truth_dir, f)
                break

    if truth_file is None:
        # Try standard location
        truth_file = os.path.join(truth_dir, 'output.sample_recession_rates.truth.csv')

    if not os.path.exists(truth_file):
        print(f"  Truth file not found in: {truth_dir}")
        return {}

    try:
        truth_df = pd.read_csv(truth_file)
        print(f"  Loaded truth from: {truth_file}")

        truth_dict = {}

        # Get column names (excluding kper)
        value_cols = [col for col in truth_df.columns if col != 'kper']

        for _, row in truth_df.iterrows():
            kper = int(row['kper']) if 'kper' in truth_df.columns else 1
            if kper not in truth_dict:
                truth_dict[kper] = {}

            for col in value_cols:
                truth_dict[kper][col] = row[col]

        return truth_dict
    except Exception as e:
        print(f"  Error loading truth file: {e}")
        return {}


def create_recession_plot(run_name, model_name, post_iter=19, filter_file=None,
                          suffix=None, pct_low=10, pct_high=90):
    """
    Create distribution plot comparing prior vs posterior for recession rates.
    """
    # Set suffix for output files
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"RECESSION RATE DISTRIBUTION ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_recession')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Get recession observations
    obs_dict, recession_df = get_recession_obs_by_kper(obs_data)

    print(f"\nRecession observations found:")
    for kper, sites in obs_dict.items():
        print(f"  kper {kper}: {list(sites.keys())}")

    # Load truth data
    print("\nLoading truth data...")
    truth_data = load_truth_data(run_name)
    if truth_data:
        for kper, sites in truth_data.items():
            print(f"  kper {kper} truth: {list(sites.keys())}")

    # Get sites (usecols) - should be same for both kper 2 and 4
    sites = []
    if 2 in obs_dict:
        sites = list(obs_dict[2].keys())
    elif 4 in obs_dict:
        sites = list(obs_dict[4].keys())

    if len(sites) == 0:
        print("Error: No recession observations found for kper 2 or 4")
        return

    print(f"\nSites: {sites}")

    # Load filtered realizations if specified
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations from {filter_file}")

    # Load prior and posterior observation ensembles
    print(f"\nLoading observation ensembles...")
    prior_file = os.path.join(master_dir, f'{model_name}.0.obs.csv')
    post_file = os.path.join(master_dir, f'{model_name}.{post_iter}.obs.csv')
    obs_noise_file = os.path.join(master_dir, f'{model_name}.obs+noise.csv')

    # Get all needed observation columns
    all_obs = []
    for kper in [2, 4]:
        if kper in obs_dict:
            all_obs.extend(obs_dict[kper].values())

    # Read only needed columns
    print(f"  Reading prior (iter 0)...")
    prior_header = pd.read_csv(prior_file, nrows=0)
    cols_to_read = [prior_header.columns[0]]  # index column
    cols_to_read.extend([col for col in all_obs if col in prior_header.columns])
    prior_oe = pd.read_csv(prior_file, usecols=cols_to_read, index_col=0)

    print(f"  Reading posterior (iter {post_iter})...")
    post_oe = pd.read_csv(post_file, usecols=cols_to_read, index_col=0)

    # Load obs+noise if available
    obs_noise_oe = None
    if os.path.exists(obs_noise_file):
        print(f"  Reading obs+noise...")
        obs_noise_header = pd.read_csv(obs_noise_file, nrows=0)
        cols_to_read_noise = [obs_noise_header.columns[0]]
        cols_to_read_noise.extend([col for col in all_obs if col in obs_noise_header.columns])
        obs_noise_oe = pd.read_csv(obs_noise_file, usecols=cols_to_read_noise, index_col=0)
        print(f"    Loaded {len(obs_noise_oe)} realizations")
    else:
        print(f"  obs+noise file not found")

    # Filter realizations if specified
    if filtered_ids:
        prior_oe = prior_oe.loc[[r for r in filtered_ids if r in prior_oe.index]]
        post_oe = post_oe.loc[[r for r in filtered_ids if r in post_oe.index]]
        if obs_noise_oe is not None:
            obs_noise_oe = obs_noise_oe.loc[[r for r in filtered_ids if r in obs_noise_oe.index]]
        print(f"  Filtered to {len(prior_oe)} prior, {len(post_oe)} posterior realizations")

    # Create figure
    print("\nCreating plot...")
    n_sites = len(sites)
    n_cols = min(n_sites, 3)
    n_rows = 2  # Row A: kper 2 (present), Row B: kper 4 (past)

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.28, 7.28/2))

    # Handle single column case
    if n_cols == 1:
        axes = axes.reshape(n_rows, 1)

    colors = {
        'prior': '#8d8d8d',      # grey
        'post_present': '#ff7f0e',  # orange
        'post_past': '#9467bd',     # purple
        'obs_noise': '#2ca02c',     # green
        'truth': 'red',
    }

    alpha = 0.6

    # Panel labels
    panel_labels_present = ['A', 'B', 'C', 'D', 'E', 'F']
    panel_labels_past = ['G', 'H', 'I', 'J', 'K', 'L']

    for site_idx, site in enumerate(sites):
        if site_idx >= n_cols:
            break

        # --- Row 0: kper 2 (present) ---
        ax_present = axes[0, site_idx]

        if 2 in obs_dict and site in obs_dict[2]:
            obs_name = obs_dict[2][site]

            if obs_name in prior_oe.columns and obs_name in post_oe.columns:
                prior_vals = prior_oe[obs_name].values
                post_vals = post_oe[obs_name].values

                # Filter to percentile range
                prior_plow, prior_phigh = np.percentile(prior_vals, [pct_low, pct_high])
                post_plow, post_phigh = np.percentile(post_vals, [pct_low, pct_high])
                prior_filt = prior_vals[(prior_vals >= prior_plow) & (prior_vals <= prior_phigh)]
                post_filt = post_vals[(post_vals >= post_plow) & (post_vals <= post_phigh)]

                # Set x-axis limits based on filtered posterior
                post_min = np.min(post_filt)
                post_max = np.max(post_filt)
                x_margin = (post_max - post_min) * 0.1
                xlim = (post_min - x_margin, post_max + x_margin)

                # Bins
                bins_prior = np.linspace(np.min(prior_filt), np.max(prior_filt), 25)
                bins_post = np.linspace(post_min, post_max, 25)

                # Plot histograms
                ax_present.hist(prior_filt, bins=bins_prior, alpha=alpha, color=colors['prior'],
                               label='Prior', edgecolor='black', linewidth=0.3)
                ax_present.hist(post_filt, bins=bins_post, alpha=alpha, color=colors['post_present'],
                               label='Posterior', edgecolor='black', linewidth=0.3)

                # Plot obs+noise if available
                median_obs_noise = None
                if obs_noise_oe is not None and obs_name in obs_noise_oe.columns:
                    obs_noise_vals = obs_noise_oe[obs_name].values
                    obs_noise_plow, obs_noise_phigh = np.percentile(obs_noise_vals, [pct_low, pct_high])
                    obs_noise_filt = obs_noise_vals[(obs_noise_vals >= obs_noise_plow) & (obs_noise_vals <= obs_noise_phigh)]
                    bins_noise = np.linspace(np.min(obs_noise_filt), np.max(obs_noise_filt), 25)
                    ax_present.hist(obs_noise_filt, bins=bins_noise, alpha=0.5, color=colors['obs_noise'],
                                   label='Obs+noise', edgecolor='black', linewidth=0.3)
                    median_obs_noise = np.median(obs_noise_vals)
                    ax_present.axvline(median_obs_noise, color=colors['obs_noise'], linestyle='-', linewidth=2)

                ax_present.set_xlim(xlim)

                # Median lines
                median_prior = np.median(prior_vals)
                median_post = np.median(post_vals)
                ax_present.axvline(median_prior, color=colors['prior'], linestyle='--', linewidth=2)
                ax_present.axvline(median_post, color=colors['post_present'], linestyle='-', linewidth=2)

                # Truth line
                truth_val = None
                if 2 in truth_data and site in truth_data[2]:
                    truth_val = truth_data[2][site]
                    ax_present.axvline(truth_val, color=colors['truth'], linestyle=':', linewidth=2, label='Truth')

                # Text box
                text_str = f"Prior: {median_prior:.4f}\nPost: {median_post:.4f}"
                if median_obs_noise is not None:
                    text_str += f"\nObs+noise: {median_obs_noise:.4f}"
                if truth_val is not None:
                    text_str += f"\nTruth: {truth_val:.4f}"
                ax_present.text(0.02, 0.98, text_str,
                               transform=ax_present.transAxes, fontsize=6, verticalalignment='top',
                               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax_present.set_xlabel('Recession rate (m/d)', fontsize=6)
        ax_present.set_ylabel('Frequency', fontsize=6)
        ax_present.set_title(f'{panel_labels_present[site_idx]}: Nov 2023 recession rate ({site})',
                            fontsize=6, fontweight='bold', loc='left')
        ax_present.grid(True, alpha=0.3)
        ax_present.tick_params(labelsize=6)

        if site_idx == 0:
            ax_present.legend(fontsize=5, loc='upper right')

        # --- Row 1: kper 4 (past) ---
        ax_past = axes[1, site_idx]

        if 4 in obs_dict and site in obs_dict[4]:
            obs_name = obs_dict[4][site]

            if obs_name in prior_oe.columns and obs_name in post_oe.columns:
                prior_vals = prior_oe[obs_name].values
                post_vals = post_oe[obs_name].values

                # Filter to percentile range
                prior_plow, prior_phigh = np.percentile(prior_vals, [pct_low, pct_high])
                post_plow, post_phigh = np.percentile(post_vals, [pct_low, pct_high])
                prior_filt = prior_vals[(prior_vals >= prior_plow) & (prior_vals <= prior_phigh)]
                post_filt = post_vals[(post_vals >= post_plow) & (post_vals <= post_phigh)]

                # Set x-axis limits based on filtered posterior
                post_min = np.min(post_filt)
                post_max = np.max(post_filt)
                x_margin = (post_max - post_min) * 0.1
                xlim = (post_min - x_margin, post_max + x_margin)

                # Bins
                bins_prior = np.linspace(np.min(prior_filt), np.max(prior_filt), 25)
                bins_post = np.linspace(post_min, post_max, 25)

                # Plot histograms
                ax_past.hist(prior_filt, bins=bins_prior, alpha=alpha, color=colors['prior'],
                            label='Prior', edgecolor='black', linewidth=0.3)
                ax_past.hist(post_filt, bins=bins_post, alpha=alpha, color=colors['post_past'],
                            label='Posterior', edgecolor='black', linewidth=0.3)

                ax_past.set_xlim(xlim)

                # Median lines
                median_prior = np.median(prior_vals)
                median_post = np.median(post_vals)
                ax_past.axvline(median_prior, color=colors['prior'], linestyle='--', linewidth=2)
                ax_past.axvline(median_post, color=colors['post_past'], linestyle='-', linewidth=2)

                # Truth line
                truth_val = None
                if 4 in truth_data and site in truth_data[4]:
                    truth_val = truth_data[4][site]
                    ax_past.axvline(truth_val, color=colors['truth'], linestyle=':', linewidth=2, label='Truth')

                # Text box
                text_str = f"Prior: {median_prior:.4f}\nPost: {median_post:.4f}"
                if truth_val is not None:
                    text_str += f"\nTruth: {truth_val:.4f}"
                ax_past.text(0.02, 0.98, text_str,
                            transform=ax_past.transAxes, fontsize=6, verticalalignment='top',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax_past.set_xlabel('Recession rate (m/d)', fontsize=6)
        ax_past.set_ylabel('Frequency', fontsize=6)
        ax_past.set_title(f'{panel_labels_past[site_idx]}: Past recession rate ({site})',
                         fontsize=6, fontweight='bold', loc='left')
        ax_past.grid(True, alpha=0.3)
        ax_past.tick_params(labelsize=6)

    # Hide unused subplots
    for site_idx in range(len(sites), n_cols):
        axes[0, site_idx].set_visible(False)
        axes[1, site_idx].set_visible(False)

    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, f'{run_name}_recession_predictions{file_suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved plot to: {output_path}")

    # Save summary statistics
    print("\nCreating summary table...")
    summary_rows = []

    for kper in [2, 4]:
        if kper not in obs_dict:
            continue

        for site, obs_name in obs_dict[kper].items():
            if obs_name not in post_oe.columns:
                continue

            prior_vals = prior_oe[obs_name].values
            post_vals = post_oe[obs_name].values

            truth_val = None
            if kper in truth_data and site in truth_data[kper]:
                truth_val = truth_data[kper][site]

            summary_rows.append({
                'kper': kper,
                'site': site,
                'obs_name': obs_name,
                'prior_min': np.min(prior_vals),
                'prior_p25': np.percentile(prior_vals, 25),
                'prior_median': np.median(prior_vals),
                'prior_mean': np.mean(prior_vals),
                'prior_p75': np.percentile(prior_vals, 75),
                'prior_max': np.max(prior_vals),
                'prior_std': np.std(prior_vals),
                'post_min': np.min(post_vals),
                'post_p25': np.percentile(post_vals, 25),
                'post_median': np.median(post_vals),
                'post_mean': np.mean(post_vals),
                'post_p75': np.percentile(post_vals, 75),
                'post_max': np.max(post_vals),
                'post_std': np.std(post_vals),
                'truth': truth_val,
            })

    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(output_dir, f'{run_name}_recession_summary{file_suffix}.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"Saved summary to: {summary_path}")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_recession_plot(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        suffix=args.suffix,
        pct_low=args.pct_low,
        pct_high=args.pct_high
    )
