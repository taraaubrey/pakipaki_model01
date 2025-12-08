"""
pp_pk4heads_predictions.py - Compare prior vs posterior pk4 head distributions

Compares prior (iter 0) and posterior distributions for pk4 head predictions
showing present vs past across different stress periods.

Supports both observation formats:
- New format: oname='pk4-heads' (from models/{run_name}/truth/output.sample_heads.truth.csv)
- Old format: oname='pk4-heads' with usecol='pk4'

Usage:
    python scripts/pp_pk4heads_predictions.py run_name [options]

Examples:
    python scripts/pp_pk4heads_predictions.py local_run34 --post-iter 19
    python scripts/pp_pk4heads_predictions.py run36 --run-name run36_62859_noGHBheads --post-iter 39
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
        description='Compare prior vs posterior pk4 head distributions'
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
    parser.add_argument('--reload', action='store_true',
                       help='Force reload data from source files (ignore cache)')
    return parser.parse_args()


def load_filtered_realizations(filter_file):
    """Load filtered realization indices from file."""
    if filter_file is None or not os.path.exists(filter_file):
        return None

    df = pd.read_csv(filter_file, comment='#')
    return df['realization'].tolist()


def get_pk4_obs_by_kper_kstp(obs_data):
    """
    Get pk4 head observation names organized by kper and kstp.

    Looks for both:
    - Old format: oname='pk4-heads' with usecol='pk4'
    - New format: oname='pk4-heads'

    Returns dict: {(kper, kstp): obs_name}
    """
    # Try new format first (oname = 'pk4-heads')
    pk4_heads = obs_data[obs_data['oname'] == 'pk4-heads'].copy()

    # If not found, try old format (pk4-heads with usecol = pk4)
    if len(pk4_heads) == 0:
        print("  No 'pk4-heads' observations found, trying 'pk4-heads' with usecol='pk4'...")
        ts_heads = obs_data[obs_data['oname'] == 'pk4-heads'].copy()
        pk4_heads = ts_heads[ts_heads['usecol'] == 'pk4'].copy()

    pk4_heads['kper'] = pd.to_numeric(pk4_heads['kper'], errors='coerce')
    pk4_heads['kstp'] = pd.to_numeric(pk4_heads['kstp'], errors='coerce')

    obs_dict = {}
    for idx, row in pk4_heads.iterrows():
        kper = int(row['kper'])
        kstp = int(row['kstp']) if pd.notna(row['kstp']) else 1
        obs_dict[(kper, kstp)] = idx

    return obs_dict, pk4_heads


def load_truth_data(truth_file):
    """
    Load truth head data for pk4.

    Time mapping:
    - time 1 = kper 1 kstp 1
    - time 2-53 = kper 2 kstp 1-52
    - time 54 = kper 3 kstp 1
    - time 55-106 = kper 4 kstp 1-52

    Returns dict: {(kper, kstp): truth_value}
    """
    if not os.path.exists(truth_file):
        print(f"  Truth file not found: {truth_file}")
        return {}

    try:
        truth_df = pd.read_csv(truth_file)

        # Check if pk4 column exists
        if 'pk4' not in truth_df.columns:
            print(f"  Warning: 'pk4' column not found in truth file")
            return {}

        truth_dict = {}

        for _, row in truth_df.iterrows():
            time_idx = int(row['time']) if 'time' in truth_df.columns else int(row.name) + 1

            # Map time to kper/kstp
            if time_idx == 1:
                kper, kstp = 1, 1
            elif 2 <= time_idx <= 53:
                kper, kstp = 2, time_idx - 1
            elif time_idx == 54:
                kper, kstp = 3, 1
            elif 55 <= time_idx <= 106:
                kper, kstp = 4, time_idx - 54
            else:
                continue

            truth_dict[(kper, kstp)] = row['pk4']

        return truth_dict
    except Exception as e:
        print(f"  Error loading truth file: {e}")
        return {}


def create_pk4_distribution_plot(run_name, model_name, post_iter=19, filter_file=None,
                                  suffix=None, pct_low=10, pct_high=90, reload=False):
    """
    Create distribution plot comparing prior vs posterior for pk4 heads.
    """
    # Set suffix for output files
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"PK4 HEADS DISTRIBUTION ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_predictions_pk4heads')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Get pk4 observations
    obs_dict, pk4_df = get_pk4_obs_by_kper_kstp(obs_data)

    print(f"\nPK4 head observations found: {len(obs_dict)}")

    # Determine available kper/kstp combinations
    kpers = sorted(set(k[0] for k in obs_dict.keys()))
    print(f"Available kper: {kpers}")

    # Get max kstp for transient periods
    kper2_kstps = sorted([k[1] for k in obs_dict.keys() if k[0] == 2])
    kper4_kstps = sorted([k[1] for k in obs_dict.keys() if k[0] == 4])

    if len(kper2_kstps) > 0:
        last_kstp_2 = max(kper2_kstps)
        print(f"kper 2: kstp 1 to {last_kstp_2}")
    else:
        last_kstp_2 = 1

    if len(kper4_kstps) > 0:
        last_kstp_4 = max(kper4_kstps)
        print(f"kper 4: kstp 1 to {last_kstp_4}")
    else:
        last_kstp_4 = 1

    # Load truth data from run-name directory
    truth_file = os.path.join('models', run_name, 'truth', 'output.sample_heads.truth.csv')
    print(f"\nLoading truth data from: {truth_file}")
    truth_data = load_truth_data(truth_file)
    print(f"  Found {len(truth_data)} truth values")

    # Define comparisons (present = kper 1-2, past = kper 3-4)
    comparisons = [
        {'name': 'Steady state', 'present': (1, 1), 'past': (3, 1)},
        {'name': '10 days after recession', 'present': (2, 10), 'past': (4, 10)},
        {'name': 'End of recession', 'present': (2, last_kstp_2), 'past': (4, last_kstp_4)},
    ]

    # Get observation names for each comparison
    for comp in comparisons:
        comp['present_obs'] = obs_dict.get(comp['present'])
        comp['past_obs'] = obs_dict.get(comp['past'])
        comp['present_truth'] = truth_data.get(comp['present'])
        comp['past_truth'] = truth_data.get(comp['past'])
        print(f"\n{comp['name']}:")
        print(f"  Present {comp['present']}: {comp['present_obs']} (truth: {comp['present_truth']})")
        print(f"  Past {comp['past']}: {comp['past_obs']} (truth: {comp['past_truth']})")

    # Load filtered realizations if specified
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations from {filter_file}")

    # Setup cache directory
    cache_dir = os.path.join(output_dir, 'cache')
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    # Cache file names
    cache_prior = os.path.join(cache_dir, f'{run_name}_prior_oe{file_suffix}.csv')
    cache_post = os.path.join(cache_dir, f'{run_name}_post_oe{file_suffix}.csv')
    cache_noise = os.path.join(cache_dir, f'{run_name}_obs_noise_oe{file_suffix}.csv')

    # Check if cached data exists
    cache_exists = os.path.exists(cache_prior) and os.path.exists(cache_post)

    # Define source file paths (needed for both cache and non-cache paths)
    prior_file = os.path.join(master_dir, f'{model_name}.0.obs.csv')
    post_file = os.path.join(master_dir, f'{model_name}.{post_iter}.obs.csv')
    obs_noise_file = os.path.join(master_dir, f'{model_name}.obs+noise.csv')

    # Read header for column information (needed for temporal data loading)
    prior_header = pd.read_csv(prior_file, nrows=0)

    if cache_exists and not reload:
        # Load from cache
        print(f"\nLoading from cache (use --reload to force refresh)...")
        print(f"  Reading cached prior...")
        prior_oe = pd.read_csv(cache_prior, index_col=0)
        print(f"  Reading cached posterior...")
        post_oe = pd.read_csv(cache_post, index_col=0)

        obs_noise_oe = None
        if os.path.exists(cache_noise):
            print(f"  Reading cached obs+noise...")
            obs_noise_oe = pd.read_csv(cache_noise, index_col=0)

        print(f"  Loaded {len(prior_oe)} prior, {len(post_oe)} posterior realizations from cache")

    else:
        # Load from source files
        print(f"\nLoading observation ensembles from source...")

        # Get all needed observation columns
        all_obs = []
        for comp in comparisons:
            if comp['present_obs']:
                all_obs.append(comp['present_obs'])
            if comp['past_obs']:
                all_obs.append(comp['past_obs'])

        # Read only needed columns
        print(f"  Reading prior (iter 0)...")
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
            print(f"  obs+noise file not found: {obs_noise_file}")

        # Filter realizations if specified
        if filtered_ids:
            prior_oe = prior_oe.loc[[r for r in filtered_ids if r in prior_oe.index]]
            post_oe = post_oe.loc[[r for r in filtered_ids if r in post_oe.index]]
            if obs_noise_oe is not None:
                obs_noise_oe = obs_noise_oe.loc[[r for r in filtered_ids if r in obs_noise_oe.index]]
            print(f"  Filtered to {len(prior_oe)} prior, {len(post_oe)} posterior realizations")

        # Save to cache
        print(f"\nSaving to cache...")
        prior_oe.to_csv(cache_prior)
        post_oe.to_csv(cache_post)
        if obs_noise_oe is not None:
            obs_noise_oe.to_csv(cache_noise)
        print(f"  Cached data saved to: {cache_dir}")

    # Save extracted data for quick reloading
    print("\nSaving extracted data...")
    data_dir = os.path.join(output_dir, 'data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Build dataframe with pk4 head data
    head_data_rows = []
    for comp in comparisons:
        present_obs = comp['present_obs']
        past_obs = comp['past_obs']
        comp_name = comp['name']

        if present_obs and present_obs in prior_oe.columns:
            for real_id in prior_oe.index:
                head_data_rows.append({
                    'comparison': comp_name,
                    'period': 'present',
                    'realization': real_id,
                    'iteration': 'prior',
                    'value': prior_oe.loc[real_id, present_obs]
                })
            for real_id in post_oe.index:
                head_data_rows.append({
                    'comparison': comp_name,
                    'period': 'present',
                    'realization': real_id,
                    'iteration': 'posterior',
                    'value': post_oe.loc[real_id, present_obs]
                })

        if past_obs and past_obs in prior_oe.columns:
            for real_id in prior_oe.index:
                head_data_rows.append({
                    'comparison': comp_name,
                    'period': 'past',
                    'realization': real_id,
                    'iteration': 'prior',
                    'value': prior_oe.loc[real_id, past_obs]
                })
            for real_id in post_oe.index:
                head_data_rows.append({
                    'comparison': comp_name,
                    'period': 'past',
                    'realization': real_id,
                    'iteration': 'posterior',
                    'value': post_oe.loc[real_id, past_obs]
                })

    head_data_df = pd.DataFrame(head_data_rows)
    head_data_file = os.path.join(data_dir, f'{run_name}_pk4heads_data{file_suffix}.csv')
    head_data_df.to_csv(head_data_file, index=False)
    print(f"  Saved: {head_data_file}")

    # Save a log of the command parameters
    import json
    log_data = {
        'run_name': run_name,
        'model_name': model_name,
        'post_iter': post_iter,
        'filter_file': filter_file,
        'suffix': suffix,
        'comparisons': [
            {'name': c['name'], 'present': c['present'], 'past': c['past']}
            for c in comparisons
        ]
    }
    log_file = os.path.join(data_dir, f'{run_name}_pk4heads_log{file_suffix}.json')
    with open(log_file, 'w') as f:
        json.dump(log_data, f, indent=2)
    print(f"  Saved: {log_file}")

    # Create figure with 3 rows: present, past, difference
    print("\nCreating plot...")
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))

    colors = {
        'prior': '#8d8d8d',      # grey
        'post_present': '#ff7f0e',  # orange (A-C)
        'post_past': '#9467bd',     # purple (D-F)
        'post_diff': '#d62728',     # red (G-I)
        'truth': '#2ca02c',         # green for truth
        'obs_noise': '#2ca02c',     # green for obs+noise
    }

    # Define titles for each panel (A-I)
    titles = {
        (0, 0): 'A: Steady state (present)',
        (0, 1): 'B: Present (10 days after recession)',
        (0, 2): 'C: Present (end of recession)',
        (1, 0): 'D: Steady state (past)',
        (1, 1): 'E: Past (10 days after recession)',
        (1, 2): 'F: Past (end of recession)',
        (2, 0): 'G: Steady state difference',
        (2, 1): 'H: 10 days difference',
        (2, 2): 'I: End of recession difference',
    }

    for col_idx, comp in enumerate(comparisons):
        present_obs = comp['present_obs']
        past_obs = comp['past_obs']
        present_truth = comp['present_truth']
        past_truth = comp['past_truth']

        # Get axes for this column
        ax_present = axes[0, col_idx]
        ax_past = axes[1, col_idx]
        ax_diff = axes[2, col_idx]

        if present_obs is None or past_obs is None:
            ax_present.text(0.5, 0.5, 'Data not available', ha='center', va='center')
            ax_past.text(0.5, 0.5, 'Data not available', ha='center', va='center')
            ax_diff.text(0.5, 0.5, 'Data not available', ha='center', va='center')
            ax_present.set_title(titles[(0, col_idx)], fontsize=6, fontweight='bold', loc='left')
            ax_past.set_title(titles[(1, col_idx)], fontsize=6, fontweight='bold', loc='left')
            ax_diff.set_title(titles[(2, col_idx)], fontsize=6, fontweight='bold', loc='left')
            continue

        # Get data
        prior_present = prior_oe[present_obs].values if present_obs in prior_oe.columns else np.array([])
        prior_past = prior_oe[past_obs].values if past_obs in prior_oe.columns else np.array([])
        post_present = post_oe[present_obs].values if present_obs in post_oe.columns else np.array([])
        post_past = post_oe[past_obs].values if past_obs in post_oe.columns else np.array([])

        # Get obs+noise data if available
        obs_noise_present = np.array([])
        if obs_noise_oe is not None and present_obs in obs_noise_oe.columns:
            obs_noise_present = obs_noise_oe[present_obs].values

        alpha = 0.6

        # --- Row 0: Present ---
        if len(post_present) > 0:
            # Filter to percentile range
            prior_plow, prior_phigh = np.percentile(prior_present, [pct_low, pct_high])
            post_plow, post_phigh = np.percentile(post_present, [pct_low, pct_high])
            prior_present_filt = prior_present[(prior_present >= prior_plow) & (prior_present <= prior_phigh)]
            post_present_filt = post_present[(post_present >= post_plow) & (post_present <= post_phigh)]

            # Set x-axis limits based on filtered posterior
            post_min = np.min(post_present_filt)
            post_max = np.max(post_present_filt)
            x_margin = (post_max - post_min) * 0.1
            xlim_present = (post_min - x_margin, post_max + x_margin)

            # Separate bins for prior and posterior
            bins_prior_present = np.linspace(np.min(prior_present_filt), np.max(prior_present_filt), 25)
            bins_post_present = np.linspace(post_min, post_max, 25)

            # Plot histograms (using filtered data)
            ax_present.hist(prior_present_filt, bins=bins_prior_present, alpha=alpha, color=colors['prior'],
                           label='Prior', edgecolor='black', linewidth=0.3)
            ax_present.hist(post_present_filt, bins=bins_post_present, alpha=alpha, color=colors['post_present'],
                           label='Posterior', edgecolor='black', linewidth=0.3)

            # Plot obs+noise histogram if available
            median_obs_noise = None
            if len(obs_noise_present) > 0:
                obs_noise_plow, obs_noise_phigh = np.percentile(obs_noise_present, [pct_low, pct_high])
                obs_noise_filt = obs_noise_present[(obs_noise_present >= obs_noise_plow) & (obs_noise_present <= obs_noise_phigh)]
                bins_obs_noise = np.linspace(np.min(obs_noise_filt), np.max(obs_noise_filt), 25)
                ax_present.hist(obs_noise_filt, bins=bins_obs_noise, alpha=0.5, color=colors['obs_noise'],
                               label='Obs+noise', edgecolor='black', linewidth=0.3)
                median_obs_noise = np.median(obs_noise_present)
                ax_present.axvline(median_obs_noise, color=colors['obs_noise'], linestyle='-', linewidth=2)

            # Set x-axis limits based on filtered posterior
            ax_present.set_xlim(xlim_present)

            # Add median lines (from full data)
            median_prior = np.median(prior_present)
            median_post = np.median(post_present)
            ax_present.axvline(median_prior, color=colors['prior'], linestyle='--', linewidth=2)
            ax_present.axvline(median_post, color=colors['post_present'], linestyle='-', linewidth=2)

            # Add truth line if available
            if present_truth is not None:
                ax_present.axvline(present_truth, color='black', linestyle=':', linewidth=2, label='Truth')

            # Add median values in top left
            text_str = f"Prior: {median_prior:.2f}\nPost: {median_post:.2f}"
            if median_obs_noise is not None:
                text_str += f"\nObs+noise: {median_obs_noise:.2f}"
            if present_truth is not None:
                text_str += f"\nTruth: {present_truth:.2f}"
            ax_present.text(0.02, 0.98, text_str,
                           transform=ax_present.transAxes, fontsize=6, verticalalignment='top',
                           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax_present.set_xlabel('Head (m)', fontsize=6)
        ax_present.set_ylabel('Frequency', fontsize=6)
        ax_present.set_title(titles[(0, col_idx)], fontsize=6, fontweight='bold', loc='left')
        ax_present.grid(True, alpha=0.3)
        ax_present.tick_params(labelsize=8)

        if col_idx == 0:
            ax_present.legend(fontsize=6, loc='upper right')

        # --- Row 1: Past ---
        if len(post_past) > 0:
            # Filter to percentile range
            prior_plow, prior_phigh = np.percentile(prior_past, [pct_low, pct_high])
            post_plow, post_phigh = np.percentile(post_past, [pct_low, pct_high])
            prior_past_filt = prior_past[(prior_past >= prior_plow) & (prior_past <= prior_phigh)]
            post_past_filt = post_past[(post_past >= post_plow) & (post_past <= post_phigh)]

            # Set x-axis limits based on filtered posterior
            post_min = np.min(post_past_filt)
            post_max = np.max(post_past_filt)
            x_margin = (post_max - post_min) * 0.1
            xlim_past = (post_min - x_margin, post_max + x_margin)

            # Separate bins for prior and posterior
            bins_prior_past = np.linspace(np.min(prior_past_filt), np.max(prior_past_filt), 25)
            bins_post_past = np.linspace(post_min, post_max, 25)

            # Plot histograms (using filtered data)
            ax_past.hist(prior_past_filt, bins=bins_prior_past, alpha=alpha, color=colors['prior'],
                        label='Prior', edgecolor='black', linewidth=0.3)
            ax_past.hist(post_past_filt, bins=bins_post_past, alpha=alpha, color=colors['post_past'],
                        label='Posterior', edgecolor='black', linewidth=0.3)

            # Set x-axis limits based on filtered posterior
            ax_past.set_xlim(xlim_past)

            # Add median lines (from full data)
            median_prior = np.median(prior_past)
            median_post = np.median(post_past)
            ax_past.axvline(median_prior, color=colors['prior'], linestyle='--', linewidth=2)
            ax_past.axvline(median_post, color=colors['post_past'], linestyle='-', linewidth=2)

            # Add truth line if available
            if past_truth is not None:
                ax_past.axvline(past_truth, color=colors['truth'], linestyle=':', linewidth=2, label='Truth')

            # Add median values in top left
            text_str = f"Prior: {median_prior:.2f}\nPost: {median_post:.2f}"
            if past_truth is not None:
                text_str += f"\nTruth: {past_truth:.2f}"
            ax_past.text(0.02, 0.98, text_str,
                        transform=ax_past.transAxes, fontsize=6, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax_past.set_xlabel('Head (m)', fontsize=6)
        ax_past.set_ylabel('Frequency', fontsize=6)
        ax_past.set_title(titles[(1, col_idx)], fontsize=6, fontweight='bold', loc='left')
        ax_past.grid(True, alpha=0.3)
        ax_past.tick_params(labelsize=8)

        # --- Row 2: Difference (past - present) ---
        if len(post_present) > 0 and len(post_past) > 0:
            # Calculate differences (past - present)
            prior_diff = prior_past - prior_present
            post_diff = post_past - post_present

            # Filter to percentile range
            prior_plow, prior_phigh = np.percentile(prior_diff, [pct_low, pct_high])
            post_plow, post_phigh = np.percentile(post_diff, [pct_low, pct_high])
            prior_diff_filt = prior_diff[(prior_diff >= prior_plow) & (prior_diff <= prior_phigh)]
            post_diff_filt = post_diff[(post_diff >= post_plow) & (post_diff <= post_phigh)]

            # Set x-axis limits based on filtered posterior difference
            post_min = np.min(post_diff_filt)
            post_max = np.max(post_diff_filt)
            x_margin = (post_max - post_min) * 0.1
            xlim_diff = (post_min - x_margin, post_max + x_margin)

            # Separate bins for prior and posterior
            bins_prior_diff = np.linspace(np.min(prior_diff_filt), np.max(prior_diff_filt), 25)
            bins_post_diff = np.linspace(post_min, post_max, 25)

            # Plot histograms (using filtered data)
            ax_diff.hist(prior_diff_filt, bins=bins_prior_diff, alpha=alpha, color=colors['prior'],
                        label='Prior', edgecolor='black', linewidth=0.3)
            ax_diff.hist(post_diff_filt, bins=bins_post_diff, alpha=alpha, color=colors['post_diff'],
                        label='Posterior', edgecolor='black', linewidth=0.3)

            # Set x-axis limits based on filtered posterior
            ax_diff.set_xlim(xlim_diff)

            # Add median lines (from full data)
            median_prior = np.median(prior_diff)
            median_post = np.median(post_diff)
            ax_diff.axvline(median_prior, color=colors['prior'], linestyle='--', linewidth=2)
            ax_diff.axvline(median_post, color=colors['post_diff'], linestyle='-', linewidth=2)

            # Add truth difference if available
            if present_truth is not None and past_truth is not None:
                truth_diff = past_truth - present_truth
                ax_diff.axvline(truth_diff, color=colors['truth'], linestyle=':', linewidth=2, label='Truth')
                text_str = f"Prior: {median_prior:.2f}\nPost: {median_post:.2f}\nTruth: {truth_diff:.2f}"
            else:
                text_str = f"Prior: {median_prior:.2f}\nPost: {median_post:.2f}"

            # Add median values in top left
            ax_diff.text(0.02, 0.98, text_str,
                        transform=ax_diff.transAxes, fontsize=6, verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax_diff.set_xlabel('Head Difference (m)', fontsize=6)
        ax_diff.set_ylabel('Frequency', fontsize=6)
        ax_diff.set_title(titles[(2, col_idx)], fontsize=6, fontweight='bold', loc='left')
        ax_diff.grid(True, alpha=0.3)
        ax_diff.tick_params(labelsize=8)

    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, f'{run_name}_pk4heads_predictions{file_suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved plot to: {output_path}")

    # --- Create temporal line plot ---
    print("\nCreating temporal line plot...")

    # Get all observations for kper 2 and kper 4
    kper2_obs = {kstp: obs_dict.get((2, kstp)) for kstp in kper2_kstps if obs_dict.get((2, kstp))}
    kper4_obs = {kstp: obs_dict.get((4, kstp)) for kstp in kper4_kstps if obs_dict.get((4, kstp))}

    if len(kper2_obs) == 0 or len(kper4_obs) == 0:
        print("  Skipping temporal plot: insufficient data")
    else:
        # Get all observation names needed
        all_temporal_obs = list(kper2_obs.values()) + list(kper4_obs.values())

        # Read additional columns if needed
        print(f"  Loading temporal data for {len(kper2_obs)} kstps (kper 2) and {len(kper4_obs)} kstps (kper 4)...")

        # Read only the temporal columns
        cols_to_read_temporal = [prior_header.columns[0]]
        cols_to_read_temporal.extend([col for col in all_temporal_obs if col in prior_header.columns])

        prior_temporal = pd.read_csv(prior_file, usecols=cols_to_read_temporal, index_col=0)
        post_temporal = pd.read_csv(post_file, usecols=cols_to_read_temporal, index_col=0)

        # Filter realizations if specified
        if filtered_ids:
            prior_temporal = prior_temporal.loc[[r for r in filtered_ids if r in prior_temporal.index]]
            post_temporal = post_temporal.loc[[r for r in filtered_ids if r in post_temporal.index]]

        # Create figure with 2 rows and 2 columns (main plot + histogram)
        from matplotlib.gridspec import GridSpec
        fig = plt.figure(figsize=(9, 7.28/2))
        gs = GridSpec(2, 2, figure=fig, width_ratios=[3, 1], wspace=0.3, hspace=0.3)

        # Panel A: kper 2 (present) - main plot
        ax_present = fig.add_subplot(gs[0, 0])
        ax_present_hist = fig.add_subplot(gs[0, 1])

        # Filter out kstp 1-2 (start from kstp 3)
        kstps_2 = sorted([k for k in kper2_obs.keys() if k >= 3])

        # Plot prior ensemble lines (grey)
        for real_id in prior_temporal.index:
            vals = [prior_temporal.loc[real_id, kper2_obs[kstp]] for kstp in kstps_2]
            ax_present.plot(kstps_2, vals, color=colors['prior'], alpha=0.15, linewidth=0.5)

        # Plot posterior ensemble lines (orange)
        for real_id in post_temporal.index:
            vals = [post_temporal.loc[real_id, kper2_obs[kstp]] for kstp in kstps_2]
            ax_present.plot(kstps_2, vals, color=colors['post_present'], alpha=0.15, linewidth=0.5)

        # Plot means (thinner, darker orange dotted for present)
        prior_mean = [np.mean(prior_temporal[kper2_obs[kstp]].values) for kstp in kstps_2]
        post_mean = [np.mean(post_temporal[kper2_obs[kstp]].values) for kstp in kstps_2]

        ax_present.plot(kstps_2, prior_mean, color='#4d4d4d', linewidth=1, linestyle=':', label='Prior mean')
        ax_present.plot(kstps_2, post_mean, color='#cc5500', linewidth=1, linestyle=':', label='Posterior mean')

        # Plot truth if available
        truth_vals_2 = [truth_data.get((2, kstp)) for kstp in kstps_2]
        if any(v is not None for v in truth_vals_2):
            valid_kstps = [k for k, v in zip(kstps_2, truth_vals_2) if v is not None]
            valid_vals = [v for v in truth_vals_2 if v is not None]
            ax_present.plot(valid_kstps, valid_vals, color='black', linewidth=2, linestyle=':', label='Truth')

        # Set y-axis limits based on truth data range with margin
        all_post_vals_2 = np.concatenate([post_temporal[kper2_obs[kstp]].values for kstp in kstps_2])
        all_prior_vals_2 = np.concatenate([prior_temporal[kper2_obs[kstp]].values for kstp in kstps_2])

        if any(v is not None for v in truth_vals_2):
            valid_truth = [v for v in truth_vals_2 if v is not None]
            truth_min = np.min(valid_truth)
            truth_max = np.max(valid_truth)
            truth_range = truth_max - truth_min
            y_min_2 = truth_min - 0.2 * truth_range
            y_max_2 = truth_max + 0.2 * truth_range
        else:
            q25_2 = np.percentile(all_post_vals_2, 25)
            q75_2 = np.percentile(all_post_vals_2, 75)
            iqr_2 = q75_2 - q25_2
            y_min_2 = max(np.min(all_post_vals_2), q25_2 - 1.5 * iqr_2)
            y_max_2 = min(np.max(all_post_vals_2), q75_2 + 1.5 * iqr_2)
        ax_present.set_ylim(y_min_2, y_max_2)

        # Create histogram for panel A (right side)
        # Collect all values across time for each ensemble
        bins_hist = np.linspace(y_min_2, y_max_2, 30)

        ax_present_hist.hist(all_prior_vals_2, bins=bins_hist, orientation='horizontal',
                            alpha=0.6, color=colors['prior'], label='Prior', edgecolor='black', linewidth=0.3)
        ax_present_hist.hist(all_post_vals_2, bins=bins_hist, orientation='horizontal',
                            alpha=0.6, color=colors['post_present'], label='Posterior', edgecolor='black', linewidth=0.3)

        # Add obs+noise histogram if available
        if obs_noise_oe is not None:
            obs_noise_vals_2 = []
            for kstp in kstps_2:
                obs_name = kper2_obs[kstp]
                if obs_name in obs_noise_oe.columns:
                    obs_noise_vals_2.extend(obs_noise_oe[obs_name].values)
            if len(obs_noise_vals_2) > 0:
                obs_noise_vals_2 = np.array(obs_noise_vals_2)
                ax_present_hist.hist(obs_noise_vals_2, bins=bins_hist, orientation='horizontal',
                                    alpha=0.5, color=colors['obs_noise'], label='Obs+noise', edgecolor='black', linewidth=0.3)

        ax_present_hist.set_ylim(y_min_2, y_max_2)
        ax_present_hist.set_xlabel('Count', fontsize=6)
        ax_present_hist.set_ylabel('')
        ax_present_hist.tick_params(labelsize=6)
        ax_present_hist.yaxis.set_ticklabels([])
        ax_present_hist.grid(True, alpha=0.3, axis='x')

        ax_present.set_xlabel('Time step (kstp)', fontsize=6)
        ax_present.set_ylabel('Head (m)', fontsize=6)
        ax_present.set_title('A: shGW1 Present (kper 2)', fontsize=6, fontweight='bold', loc='left')
        ax_present.grid(True, alpha=0.3)
        ax_present.legend(fontsize=6, loc='best')
        ax_present.tick_params(labelsize=6)

        # Panel B: kper 4 (past) - main plot
        ax_past = fig.add_subplot(gs[1, 0])
        ax_past_hist = fig.add_subplot(gs[1, 1])

        # Filter out kstp 1-2 (start from kstp 3)
        kstps_4 = sorted([k for k in kper4_obs.keys() if k >= 3])

        # Plot prior ensemble lines (grey)
        for real_id in prior_temporal.index:
            vals = [prior_temporal.loc[real_id, kper4_obs[kstp]] for kstp in kstps_4]
            ax_past.plot(kstps_4, vals, color=colors['prior'], alpha=0.15, linewidth=0.5)

        # Plot posterior ensemble lines (purple)
        for real_id in post_temporal.index:
            vals = [post_temporal.loc[real_id, kper4_obs[kstp]] for kstp in kstps_4]
            ax_past.plot(kstps_4, vals, color=colors['post_past'], alpha=0.15, linewidth=0.5)

        # Plot means (thinner, darker purple dotted for past)
        prior_mean_4 = [np.mean(prior_temporal[kper4_obs[kstp]].values) for kstp in kstps_4]
        post_mean_4 = [np.mean(post_temporal[kper4_obs[kstp]].values) for kstp in kstps_4]

        ax_past.plot(kstps_4, prior_mean_4, color='#4d4d4d', linewidth=1, linestyle=':', label='Prior mean')
        ax_past.plot(kstps_4, post_mean_4, color='#5c3d6e', linewidth=1, linestyle=':', label='Posterior mean')

        # Plot truth if available
        truth_vals_4 = [truth_data.get((4, kstp)) for kstp in kstps_4]
        if any(v is not None for v in truth_vals_4):
            valid_kstps = [k for k, v in zip(kstps_4, truth_vals_4) if v is not None]
            valid_vals = [v for v in truth_vals_4 if v is not None]
            ax_past.plot(valid_kstps, valid_vals, color='black', linewidth=2, linestyle=':', label='Truth')

        # Set y-axis limits based on IQR * 1.5 of posterior values
        all_post_vals_4 = np.concatenate([post_temporal[kper4_obs[kstp]].values for kstp in kstps_4])
        all_prior_vals_4 = np.concatenate([prior_temporal[kper4_obs[kstp]].values for kstp in kstps_4])

        q25_4 = np.percentile(all_post_vals_4, 25)
        q75_4 = np.percentile(all_post_vals_4, 75)
        iqr_4 = q75_4 - q25_4
        y_min_4 = max(np.min(all_post_vals_4), q25_4 - 1.5 * iqr_4)
        y_max_4 = min(np.max(all_post_vals_4), q75_4 + 1.5 * iqr_4)
        ax_past.set_ylim(y_min_4, y_max_4)

        # Create histogram for panel B (right side)
        bins_hist_4 = np.linspace(y_min_4, y_max_4, 30)

        ax_past_hist.hist(all_prior_vals_4, bins=bins_hist_4, orientation='horizontal',
                         alpha=0.6, color=colors['prior'], label='Prior', edgecolor='black', linewidth=0.3)
        ax_past_hist.hist(all_post_vals_4, bins=bins_hist_4, orientation='horizontal',
                         alpha=0.6, color=colors['post_past'], label='Posterior', edgecolor='black', linewidth=0.3)

        # Add obs+noise histogram if available
        if obs_noise_oe is not None:
            obs_noise_vals_4 = []
            for kstp in kstps_4:
                obs_name = kper4_obs[kstp]
                if obs_name in obs_noise_oe.columns:
                    obs_noise_vals_4.extend(obs_noise_oe[obs_name].values)
            if len(obs_noise_vals_4) > 0:
                obs_noise_vals_4 = np.array(obs_noise_vals_4)
                ax_past_hist.hist(obs_noise_vals_4, bins=bins_hist_4, orientation='horizontal',
                                 alpha=0.5, color=colors['obs_noise'], label='Obs+noise', edgecolor='black', linewidth=0.3)

        ax_past_hist.set_ylim(y_min_4, y_max_4)
        ax_past_hist.set_xlabel('Count', fontsize=6)
        ax_past_hist.set_ylabel('')
        ax_past_hist.tick_params(labelsize=6)
        ax_past_hist.yaxis.set_ticklabels([])
        ax_past_hist.grid(True, alpha=0.3, axis='x')

        ax_past.set_xlabel('Time step (kstp)', fontsize=6)
        ax_past.set_ylabel('Head (m)', fontsize=6)
        ax_past.set_title('B: shGW1 Past (kper 4)', fontsize=6, fontweight='bold', loc='left')
        ax_past.grid(True, alpha=0.3)
        ax_past.legend(fontsize=6, loc='best')
        ax_past.tick_params(labelsize=6)

        plt.tight_layout()

        # Save temporal plot
        temporal_path = os.path.join(output_dir, f'{run_name}_shGW1_temporal{file_suffix}.png')
        plt.savefig(temporal_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"Saved temporal plot to: {temporal_path}")

        # --- Create summary table ---
        print("\nCreating summary table...")

        summary_rows = []

        # Process kper 2 (present)
        for kstp in kstps_2:
            if kstp == 1 or (kstp - 1) % 10 == 0:  # 1, 11, 21, 31, ...
                obs_name = kper2_obs[kstp]
                prior_vals = prior_temporal[obs_name].values
                post_vals = post_temporal[obs_name].values
                truth_val = truth_data.get((2, kstp))

                summary_rows.append({
                    'kper': 2,
                    'kstp': kstp,
                    'truth': truth_val,
                    'prior_mean': np.mean(prior_vals),
                    'prior_std': np.std(prior_vals),
                    'prior_min': np.min(prior_vals),
                    'prior_p10': np.percentile(prior_vals, 10),
                    'prior_p50': np.percentile(prior_vals, 50),
                    'prior_p90': np.percentile(prior_vals, 90),
                    'prior_max': np.max(prior_vals),
                    'post_mean': np.mean(post_vals),
                    'post_std': np.std(post_vals),
                    'post_min': np.min(post_vals),
                    'post_p10': np.percentile(post_vals, 10),
                    'post_p50': np.percentile(post_vals, 50),
                    'post_p90': np.percentile(post_vals, 90),
                    'post_max': np.max(post_vals),
                })

        # Process kper 4 (past)
        for kstp in kstps_4:
            if kstp == 1 or (kstp - 1) % 10 == 0:  # 1, 11, 21, 31, ...
                obs_name = kper4_obs[kstp]
                prior_vals = prior_temporal[obs_name].values
                post_vals = post_temporal[obs_name].values
                truth_val = truth_data.get((4, kstp))

                summary_rows.append({
                    'kper': 4,
                    'kstp': kstp,
                    'truth': truth_val,
                    'prior_mean': np.mean(prior_vals),
                    'prior_std': np.std(prior_vals),
                    'prior_min': np.min(prior_vals),
                    'prior_p10': np.percentile(prior_vals, 10),
                    'prior_p50': np.percentile(prior_vals, 50),
                    'prior_p90': np.percentile(prior_vals, 90),
                    'prior_max': np.max(prior_vals),
                    'post_mean': np.mean(post_vals),
                    'post_std': np.std(post_vals),
                    'post_min': np.min(post_vals),
                    'post_p10': np.percentile(post_vals, 10),
                    'post_p50': np.percentile(post_vals, 50),
                    'post_p90': np.percentile(post_vals, 90),
                    'post_max': np.max(post_vals),
                })

        summary_df = pd.DataFrame(summary_rows)
        summary_file = os.path.join(data_dir, f'{run_name}_pk4heads_summary{file_suffix}.csv')
        summary_df.to_csv(summary_file, index=False)
        print(f"Saved summary table to: {summary_file}")

        # --- Calculate head decline statistics ---
        print("\nCalculating head decline statistics...")

        decline_stats = []

        # kper 2: decline from kstp 1 to kstp 10
        if 1 in kper2_obs and 10 in kper2_obs:
            obs_kstp1 = kper2_obs[1]
            obs_kstp10 = kper2_obs[10]

            if obs_kstp1 in prior_temporal.columns and obs_kstp10 in prior_temporal.columns:
                # Decline = kstp1 - kstp10 (positive means head dropped)
                prior_decline = prior_temporal[obs_kstp1].values - prior_temporal[obs_kstp10].values
                post_decline = post_temporal[obs_kstp1].values - post_temporal[obs_kstp10].values

                decline_stats.append({
                    'kper': 2,
                    'from_kstp': 1,
                    'to_kstp': 10,
                    'description': 'Present (kper 2) decline kstp 1 to 10',
                    'prior_min': np.min(prior_decline),
                    'prior_p05': np.percentile(prior_decline, 5),
                    'prior_p25': np.percentile(prior_decline, 25),
                    'prior_median': np.median(prior_decline),
                    'prior_mean': np.mean(prior_decline),
                    'prior_p75': np.percentile(prior_decline, 75),
                    'prior_p95': np.percentile(prior_decline, 95),
                    'prior_max': np.max(prior_decline),
                    'prior_std': np.std(prior_decline),
                    'post_min': np.min(post_decline),
                    'post_p05': np.percentile(post_decline, 5),
                    'post_p25': np.percentile(post_decline, 25),
                    'post_median': np.median(post_decline),
                    'post_mean': np.mean(post_decline),
                    'post_p75': np.percentile(post_decline, 75),
                    'post_p95': np.percentile(post_decline, 95),
                    'post_max': np.max(post_decline),
                    'post_std': np.std(post_decline),
                })

                print(f"\n  kper 2 (Present) head decline (kstp 1 to 10):")
                print(f"    Prior:     median={np.median(prior_decline):.3f}, mean={np.mean(prior_decline):.3f}, std={np.std(prior_decline):.3f}")
                print(f"    Posterior: median={np.median(post_decline):.3f}, mean={np.mean(post_decline):.3f}, std={np.std(post_decline):.3f}")
                print(f"    Posterior range: [{np.min(post_decline):.3f}, {np.max(post_decline):.3f}]")
                print(f"    Posterior IQR: [{np.percentile(post_decline, 25):.3f}, {np.percentile(post_decline, 75):.3f}]")

        # kper 4: decline from kstp 1 to kstp 10
        if 1 in kper4_obs and 10 in kper4_obs:
            obs_kstp1 = kper4_obs[1]
            obs_kstp10 = kper4_obs[10]

            if obs_kstp1 in prior_temporal.columns and obs_kstp10 in prior_temporal.columns:
                prior_decline = prior_temporal[obs_kstp1].values - prior_temporal[obs_kstp10].values
                post_decline = post_temporal[obs_kstp1].values - post_temporal[obs_kstp10].values

                decline_stats.append({
                    'kper': 4,
                    'from_kstp': 1,
                    'to_kstp': 10,
                    'description': 'Past (kper 4) decline kstp 1 to 10',
                    'prior_min': np.min(prior_decline),
                    'prior_p05': np.percentile(prior_decline, 5),
                    'prior_p25': np.percentile(prior_decline, 25),
                    'prior_median': np.median(prior_decline),
                    'prior_mean': np.mean(prior_decline),
                    'prior_p75': np.percentile(prior_decline, 75),
                    'prior_p95': np.percentile(prior_decline, 95),
                    'prior_max': np.max(prior_decline),
                    'prior_std': np.std(prior_decline),
                    'post_min': np.min(post_decline),
                    'post_p05': np.percentile(post_decline, 5),
                    'post_p25': np.percentile(post_decline, 25),
                    'post_median': np.median(post_decline),
                    'post_mean': np.mean(post_decline),
                    'post_p75': np.percentile(post_decline, 75),
                    'post_p95': np.percentile(post_decline, 95),
                    'post_max': np.max(post_decline),
                    'post_std': np.std(post_decline),
                })

                print(f"\n  kper 4 (Past) head decline (kstp 1 to 10):")
                print(f"    Prior:     median={np.median(prior_decline):.3f}, mean={np.mean(prior_decline):.3f}, std={np.std(prior_decline):.3f}")
                print(f"    Posterior: median={np.median(post_decline):.3f}, mean={np.mean(post_decline):.3f}, std={np.std(post_decline):.3f}")
                print(f"    Posterior range: [{np.min(post_decline):.3f}, {np.max(post_decline):.3f}]")
                print(f"    Posterior IQR: [{np.percentile(post_decline, 25):.3f}, {np.percentile(post_decline, 75):.3f}]")

        if decline_stats:
            decline_df = pd.DataFrame(decline_stats)
            decline_file = os.path.join(data_dir, f'{run_name}_pk4heads_decline{file_suffix}.csv')
            decline_df.to_csv(decline_file, index=False)
            print(f"\nSaved head decline statistics to: {decline_file}")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_pk4_distribution_plot(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        suffix=args.suffix,
        pct_low=args.pct_low,
        pct_high=args.pct_high,
        reload=args.reload
    )
