"""
q_spring_flux_distributions.py - Compare prior vs posterior spring flux distributions

Compares prior (iter 0) and posterior distributions for spring flux predictions
showing present vs future across different stress periods.

Usage:
    python scripts/q_spring_flux_distributions.py run_name [options]

Examples:
    python scripts/q_spring_flux_distributions.py local_run34 --post-iter 19
    python scripts/o_filter_realizations.py local_run34 --run-name 62730_run34_rch250_rr -i 19 --max-head 20 --no-positive-confined-flux
    python scripts/q_spring_flux_distributions.py local_run34 --post-iter 19 --filter-file models/local_run34/filtered_realizations.csv
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
        description='Compare prior vs posterior spring flux distributions'
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
    parser.add_argument('--spring-i', type=int, default=12,
                       help='Spring cell i index (default: 12)')
    parser.add_argument('--spring-j', type=int, default=18,
                       help='Spring cell j index (default: 18)')
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


def get_spring_obs_by_kper_kstp(obs_data, spring_i, spring_j):
    """
    Get spring flux observation names organized by kper and kstp.

    Supports both formats:
    - Old: oname='arr-spq' (array observations)
    - New: oname='spring-q' (from GHB_SP_fluxes)

    Returns dict: {(kper, kstp): obs_name}
    """
    # Try new format first (oname = 'spring-q')
    spring_q = obs_data[obs_data['oname'] == 'spring-q'].copy()

    if len(spring_q) > 0:
        spring_q['i'] = pd.to_numeric(spring_q['i'], errors='coerce')
        spring_q['j'] = pd.to_numeric(spring_q['j'], errors='coerce')
        spring_q['kper'] = pd.to_numeric(spring_q['kper'], errors='coerce')
        spring_q['kstp'] = pd.to_numeric(spring_q['kstp'], errors='coerce')

        spring = spring_q[(spring_q['i'] == spring_i) & (spring_q['j'] == spring_j)]
    else:
        # Fall back to old format (oname = 'arr-spq')
        print("  No 'spring-q' observations found, trying 'arr-spq'...")
        arr_spq = obs_data[obs_data['oname'] == 'arr-spq'].copy()
        arr_spq['i'] = pd.to_numeric(arr_spq['i'], errors='coerce')
        arr_spq['j'] = pd.to_numeric(arr_spq['j'], errors='coerce')
        arr_spq['kper'] = pd.to_numeric(arr_spq['kper'], errors='coerce')
        arr_spq['kstp'] = pd.to_numeric(arr_spq['kstp'], errors='coerce')

        spring = arr_spq[(arr_spq['i'] == spring_i) & (arr_spq['j'] == spring_j)]

    obs_dict = {}
    for idx, row in spring.iterrows():
        kper = int(row['kper'])
        kstp = int(row['kstp']) if pd.notna(row['kstp']) else 1
        obs_dict[(kper, kstp)] = idx

    return obs_dict, spring


def create_spring_distribution_plot(run_name, model_name, post_iter=19, filter_file=None,
                                    spring_i=12, spring_j=18, suffix=None, pct_low=10, pct_high=90):
    """
    Create distribution plot comparing prior vs posterior for spring flux.
    """
    # Set suffix for output files
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"SPRING FLUX DISTRIBUTION ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_predictions_springflux')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Get spring observations
    obs_dict, spring_df = get_spring_obs_by_kper_kstp(obs_data, spring_i, spring_j)

    print(f"\nSpring flux observations at i={spring_i}, j={spring_j}")

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
        print(f"\n{comp['name']}:")
        print(f"  Present {comp['present']}: {comp['present_obs']}")
        print(f"  Past {comp['past']}: {comp['past_obs']}")

    # Load filtered realizations if specified
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations from {filter_file}")

    # Load prior and posterior observation ensembles
    print(f"\nLoading observation ensembles...")
    prior_file = os.path.join(master_dir, f'{model_name}.0.obs.csv')
    post_file = os.path.join(master_dir, f'{model_name}.{post_iter}.obs.csv')

    # Get all needed observation columns
    all_obs = []
    for comp in comparisons:
        if comp['present_obs']:
            all_obs.append(comp['present_obs'])
        if comp['past_obs']:
            all_obs.append(comp['past_obs'])

    # Read only needed columns
    print(f"  Reading prior (iter 0)...")
    prior_header = pd.read_csv(prior_file, nrows=0)
    cols_to_read = [prior_header.columns[0]]  # index column
    cols_to_read.extend([col for col in all_obs if col in prior_header.columns])
    prior_oe = pd.read_csv(prior_file, usecols=cols_to_read, index_col=0)

    print(f"  Reading posterior (iter {post_iter})...")
    post_oe = pd.read_csv(post_file, usecols=cols_to_read, index_col=0)

    # Filter realizations if specified
    if filtered_ids:
        prior_oe = prior_oe.loc[[r for r in filtered_ids if r in prior_oe.index]]
        post_oe = post_oe.loc[[r for r in filtered_ids if r in post_oe.index]]
        print(f"  Filtered to {len(prior_oe)} prior, {len(post_oe)} posterior realizations")

    # Save extracted data for quick reloading
    print("\nSaving extracted data...")
    data_dir = os.path.join(output_dir, 'data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Build dataframe with spring flux data
    spring_data_rows = []
    for comp in comparisons:
        present_obs = comp['present_obs']
        past_obs = comp['past_obs']
        comp_name = comp['name']

        if present_obs and present_obs in prior_oe.columns:
            for real_id in prior_oe.index:
                spring_data_rows.append({
                    'comparison': comp_name,
                    'period': 'present',
                    'realization': real_id,
                    'iteration': 'prior',
                    'value': prior_oe.loc[real_id, present_obs]
                })
            for real_id in post_oe.index:
                spring_data_rows.append({
                    'comparison': comp_name,
                    'period': 'present',
                    'realization': real_id,
                    'iteration': 'posterior',
                    'value': post_oe.loc[real_id, present_obs]
                })

        if past_obs and past_obs in prior_oe.columns:
            for real_id in prior_oe.index:
                spring_data_rows.append({
                    'comparison': comp_name,
                    'period': 'past',
                    'realization': real_id,
                    'iteration': 'prior',
                    'value': prior_oe.loc[real_id, past_obs]
                })
            for real_id in post_oe.index:
                spring_data_rows.append({
                    'comparison': comp_name,
                    'period': 'past',
                    'realization': real_id,
                    'iteration': 'posterior',
                    'value': post_oe.loc[real_id, past_obs]
                })

    spring_data_df = pd.DataFrame(spring_data_rows)
    spring_data_file = os.path.join(data_dir, f'{run_name}_spring_flux_data{file_suffix}.csv')
    spring_data_df.to_csv(spring_data_file, index=False)
    print(f"  Saved: {spring_data_file}")

    # Save a log of the command parameters
    import json
    log_data = {
        'run_name': run_name,
        'model_name': model_name,
        'post_iter': post_iter,
        'filter_file': filter_file,
        'spring_i': spring_i,
        'spring_j': spring_j,
        'suffix': suffix,
        'comparisons': [
            {'name': c['name'], 'present': c['present'], 'past': c['past']}
            for c in comparisons
        ]
    }
    log_file = os.path.join(data_dir, f'{run_name}_spring_flux_log{file_suffix}.json')
    with open(log_file, 'w') as f:
        json.dump(log_data, f, indent=2)
    print(f"  Saved: {log_file}")

    # Create figure with 3 rows: present, past, difference
    print("\nCreating plot...")
    fig, axes = plt.subplots(3, 3, figsize=(7.28, 7.28))

    colors = {
        'prior': '#8d8d8d',      # grey
        'post_present': '#ff7f0e',  # orange (A-C)
        'post_past': '#9467bd',     # purple (D-F)
        'post_diff': '#d62728',     # red (G-I)
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

            # Set x-axis limits based on filtered posterior
            ax_present.set_xlim(xlim_present)

        ax_present.set_xlabel('Spring Flux (m³/d)', fontsize=9)
        ax_present.set_ylabel('Frequency', fontsize=9)
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

        ax_past.set_xlabel('Spring Flux (m³/d)', fontsize=9)
        ax_past.set_ylabel('Frequency', fontsize=9)
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

        ax_diff.set_xlabel('Flux Difference (m³/d)', fontsize=9)
        ax_diff.set_ylabel('Frequency', fontsize=9)
        ax_diff.set_title(titles[(2, col_idx)], fontsize=6, fontweight='bold', loc='left')
        ax_diff.grid(True, alpha=0.3)
        ax_diff.tick_params(labelsize=8)

    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, f'{run_name}_spring_flux_predictions{file_suffix}.png')
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

        # Create figure with 2 rows
        fig, axes = plt.subplots(2, 1, figsize=(7.28, 7.28/2))

        # Panel A: kper 2 (present)
        ax_present = axes[0]
        # Filter out kstp 1-2 (start from kstp 3)
        kstps_2 = sorted([k for k in kper2_obs.keys() if k >= 3])

        # Plot prior ensemble lines (grey)
        for real_id in prior_temporal.index:
            vals = [prior_temporal.loc[real_id, kper2_obs[kstp]] for kstp in kstps_2]
            ax_present.plot(kstps_2, vals, color=colors['prior'], alpha=0.15, linewidth=0.5)

        # Plot posterior ensemble lines (orange)
        for real_id in post_temporal.index:
            vals = [post_temporal.loc[real_id, kper2_obs[kstp]] for kstp in kstps_2]
            ax_present.plot(kstps_2, vals, color=colors['post_present'], alpha=0.4, linewidth=0.5)

        # Calculate statistics for posterior
        post_p10 = [np.percentile(post_temporal[kper2_obs[kstp]].values, 10) for kstp in kstps_2]
        post_p50 = [np.median(post_temporal[kper2_obs[kstp]].values) for kstp in kstps_2]
        post_p90 = [np.percentile(post_temporal[kper2_obs[kstp]].values, 90) for kstp in kstps_2]
        post_mean = [np.mean(post_temporal[kper2_obs[kstp]].values) for kstp in kstps_2]
        post_std = [np.std(post_temporal[kper2_obs[kstp]].values) for kstp in kstps_2]
        post_mean_plus_std = [m + s for m, s in zip(post_mean, post_std)]
        post_mean_minus_std = [m - s for m, s in zip(post_mean, post_std)]

        # Calculate statistics for prior
        prior_p50 = [np.median(prior_temporal[kper2_obs[kstp]].values) for kstp in kstps_2]

        # Plot prior median
        ax_present.plot(kstps_2, prior_p50, color=colors['prior'], linewidth=1.5, linestyle='--', label='Prior median')

        # Plot posterior statistics
        ax_present.plot(kstps_2, post_p50, color=colors['post_present'], linewidth=2, label='Posterior median')
        ax_present.plot(kstps_2, post_p10, color=colors['post_present'], linewidth=1, linestyle=':', label='P10/P90')
        ax_present.plot(kstps_2, post_p90, color=colors['post_present'], linewidth=1, linestyle=':')
        ax_present.plot(kstps_2, post_mean_plus_std, color=colors['post_present'], linewidth=1, linestyle='-.', alpha=0.7, label='Mean±1σ')
        ax_present.plot(kstps_2, post_mean_minus_std, color=colors['post_present'], linewidth=1, linestyle='-.', alpha=0.7)

        # Set y-axis limits based on posterior values
        all_post_vals_2 = np.concatenate([post_temporal[kper2_obs[kstp]].values for kstp in kstps_2])
        ax_present.set_ylim(np.min(all_post_vals_2), np.max(all_post_vals_2))

        ax_present.set_xlabel('Time step (kstp)', fontsize=6)
        ax_present.set_ylabel('Spring Flux (m³/d)', fontsize=6)
        ax_present.set_title('A: Present (kper 2)', fontsize=6, fontweight='bold', loc='left')
        ax_present.grid(True, alpha=0.3)
        ax_present.legend(fontsize=6, loc='best')
        ax_present.tick_params(labelsize=8)

        # Panel B: kper 4 (past)
        ax_past = axes[1]
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

        # Calculate statistics for posterior
        post_p10_4 = [np.percentile(post_temporal[kper4_obs[kstp]].values, 10) for kstp in kstps_4]
        post_p50_4 = [np.median(post_temporal[kper4_obs[kstp]].values) for kstp in kstps_4]
        post_p90_4 = [np.percentile(post_temporal[kper4_obs[kstp]].values, 90) for kstp in kstps_4]
        post_mean_4 = [np.mean(post_temporal[kper4_obs[kstp]].values) for kstp in kstps_4]
        post_std_4 = [np.std(post_temporal[kper4_obs[kstp]].values) for kstp in kstps_4]
        post_mean_plus_std_4 = [m + s for m, s in zip(post_mean_4, post_std_4)]
        post_mean_minus_std_4 = [m - s for m, s in zip(post_mean_4, post_std_4)]

        # Calculate statistics for prior
        prior_p50_4 = [np.median(prior_temporal[kper4_obs[kstp]].values) for kstp in kstps_4]

        # Plot prior median
        ax_past.plot(kstps_4, prior_p50_4, color=colors['prior'], linewidth=1.5, linestyle='--', label='Prior median')

        # Plot posterior statistics
        ax_past.plot(kstps_4, post_p50_4, color=colors['post_past'], linewidth=2, label='Posterior median')
        ax_past.plot(kstps_4, post_p10_4, color=colors['post_past'], linewidth=1, linestyle=':', label='P10/P90')
        ax_past.plot(kstps_4, post_p90_4, color=colors['post_past'], linewidth=1, linestyle=':')
        ax_past.plot(kstps_4, post_mean_plus_std_4, color=colors['post_past'], linewidth=1, linestyle='-.', alpha=0.7, label='Mean±1σ')
        ax_past.plot(kstps_4, post_mean_minus_std_4, color=colors['post_past'], linewidth=1, linestyle='-.', alpha=0.7)

        # Set y-axis limits based on posterior values
        all_post_vals_4 = np.concatenate([post_temporal[kper4_obs[kstp]].values for kstp in kstps_4])
        ax_past.set_ylim(np.min(all_post_vals_4), np.max(all_post_vals_4))

        ax_past.set_xlabel('Time step (kstp)', fontsize=6)
        ax_past.set_ylabel('Spring Flux (m³/d)', fontsize=6)
        ax_past.set_title('B: Past (kper 4)', fontsize=6, fontweight='bold', loc='left')
        ax_past.grid(True, alpha=0.3)
        ax_past.legend(fontsize=6, loc='best')
        ax_past.tick_params(labelsize=8)

        plt.tight_layout()

        # Save temporal plot
        temporal_path = os.path.join(output_dir, f'{run_name}_spring_flux_temporal{file_suffix}.png')
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

                summary_rows.append({
                    'kper': 2,
                    'kstp': kstp,
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

                summary_rows.append({
                    'kper': 4,
                    'kstp': kstp,
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
        summary_file = os.path.join(data_dir, f'{run_name}_spring_flux_summary{file_suffix}.csv')
        summary_df.to_csv(summary_file, index=False)
        print(f"Saved summary table to: {summary_file}")

    # --- Create key timesteps statistics table ---
    print("\nCreating key timesteps statistics table...")

    # Define key timesteps
    key_timesteps = [
        (1, 1),
        (2, 10),
        (2, 52),
        (3, 1),
        (4, 10),
        (4, 52),
    ]

    # Need to load data for all these timesteps
    key_obs_names = []
    for kper, kstp in key_timesteps:
        obs_name = obs_dict.get((kper, kstp))
        if obs_name:
            key_obs_names.append(obs_name)

    if len(key_obs_names) > 0:
        # Read the required columns
        cols_to_read_key = [prior_header.columns[0]]
        cols_to_read_key.extend([col for col in key_obs_names if col in prior_header.columns])

        prior_key = pd.read_csv(prior_file, usecols=cols_to_read_key, index_col=0)
        post_key = pd.read_csv(post_file, usecols=cols_to_read_key, index_col=0)

        # Filter realizations if specified
        if filtered_ids:
            prior_key = prior_key.loc[[r for r in filtered_ids if r in prior_key.index]]
            post_key = post_key.loc[[r for r in filtered_ids if r in post_key.index]]

        # Build statistics table
        key_stats_rows = []
        for kper, kstp in key_timesteps:
            obs_name = obs_dict.get((kper, kstp))
            if obs_name and obs_name in post_key.columns:
                prior_vals = prior_key[obs_name].values
                post_vals = post_key[obs_name].values

                key_stats_rows.append({
                    'kper': kper,
                    'kstp': kstp,
                    'prior_min': np.min(prior_vals),
                    'prior_max': np.max(prior_vals),
                    'prior_mean': np.mean(prior_vals),
                    'prior_median': np.median(prior_vals),
                    'post_min': np.min(post_vals),
                    'post_max': np.max(post_vals),
                    'post_mean': np.mean(post_vals),
                    'post_median': np.median(post_vals),
                })
            else:
                key_stats_rows.append({
                    'kper': kper,
                    'kstp': kstp,
                    'prior_min': np.nan,
                    'prior_max': np.nan,
                    'prior_mean': np.nan,
                    'prior_median': np.nan,
                    'post_min': np.nan,
                    'post_max': np.nan,
                    'post_mean': np.nan,
                    'post_median': np.nan,
                })

        key_stats_df = pd.DataFrame(key_stats_rows)
        key_stats_file = os.path.join(output_dir, f'{run_name}_spring_flux_key_stats{file_suffix}.csv')
        key_stats_df.to_csv(key_stats_file, index=False)
        print(f"Saved key timesteps statistics to: {key_stats_file}")
    else:
        print("  No key timesteps data available")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_spring_distribution_plot(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        spring_i=args.spring_i,
        spring_j=args.spring_j,
        suffix=args.suffix,
        pct_low=args.pct_low,
        pct_high=args.pct_high
    )
