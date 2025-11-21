"""
q_spring_flux_distributions.py - Compare prior vs posterior spring flux distributions

Compares prior (iter 0) and posterior distributions for spring flux predictions
showing present vs future across different stress periods.

Usage:
    python scripts/q_spring_flux_distributions.py run_name [options]

Examples:
    python scripts/q_spring_flux_distributions.py local_run34 --post-iter 19
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


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Compare prior vs posterior spring flux distributions'
    )
    parser.add_argument('run_name', type=str, help='Name of the run directory')
    parser.add_argument('--model-name', type=str, default=None,
                       help='Model name (default: same as run_name)')
    parser.add_argument('--post-iter', type=int, default=19,
                       help='Posterior iteration (default: 19)')
    parser.add_argument('--filter-file', type=str, default=None,
                       help='Path to file with filtered realization indices')
    parser.add_argument('--spring-i', type=int, default=36,
                       help='Spring cell i index (default: 36)')
    parser.add_argument('--spring-j', type=int, default=33,
                       help='Spring cell j index (default: 33)')
    parser.add_argument('--suffix', type=str, default=None,
                       help='Suffix for output files (e.g., "_filtered")')
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

    Returns dict: {(kper, kstp): obs_name}
    """
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
                                    spring_i=36, spring_j=33, suffix=None):
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
    output_dir = os.path.join('models', run_name, 'pp_prediction_springflux')

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

    # Define comparisons
    comparisons = [
        {'name': 'A. Steady state', 'present': (1, 1), 'future': (3, 1)},
        {'name': 'B. Transient (first kstp)', 'present': (2, 1), 'future': (4, 1)},
        {'name': f'C. Transient (last kstp)', 'present': (2, last_kstp_2), 'future': (4, last_kstp_4)},
    ]

    # Get observation names for each comparison
    for comp in comparisons:
        comp['present_obs'] = obs_dict.get(comp['present'])
        comp['future_obs'] = obs_dict.get(comp['future'])
        print(f"\n{comp['name']}:")
        print(f"  Present {comp['present']}: {comp['present_obs']}")
        print(f"  Future {comp['future']}: {comp['future_obs']}")

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
        if comp['future_obs']:
            all_obs.append(comp['future_obs'])

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
        future_obs = comp['future_obs']
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

        if future_obs and future_obs in prior_oe.columns:
            for real_id in prior_oe.index:
                spring_data_rows.append({
                    'comparison': comp_name,
                    'period': 'future',
                    'realization': real_id,
                    'iteration': 'prior',
                    'value': prior_oe.loc[real_id, future_obs]
                })
            for real_id in post_oe.index:
                spring_data_rows.append({
                    'comparison': comp_name,
                    'period': 'future',
                    'realization': real_id,
                    'iteration': 'posterior',
                    'value': post_oe.loc[real_id, future_obs]
                })

    spring_data_df = pd.DataFrame(spring_data_rows)
    spring_data_file = os.path.join(data_dir, f'{run_name}_spring_flux_data{file_suffix}.csv')
    spring_data_df.to_csv(spring_data_file, index=False)
    print(f"  Saved: {spring_data_file}")

    # Create figure
    print("\nCreating plot...")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    colors = {
        'prior_present': '#1f77b4',    # blue
        'prior_future': '#aec7e8',     # light blue
        'post_present': '#ff7f0e',     # orange
        'post_future': '#ffbb78',      # light orange
    }

    for idx, (ax, comp) in enumerate(zip(axes, comparisons)):
        present_obs = comp['present_obs']
        future_obs = comp['future_obs']

        if present_obs is None or future_obs is None:
            ax.text(0.5, 0.5, 'Data not available', ha='center', va='center')
            ax.set_title(comp['name'])
            continue

        # Get data
        prior_present = prior_oe[present_obs].values if present_obs in prior_oe.columns else []
        prior_future = prior_oe[future_obs].values if future_obs in prior_oe.columns else []
        post_present = post_oe[present_obs].values if present_obs in post_oe.columns else []
        post_future = post_oe[future_obs].values if future_obs in post_oe.columns else []

        # Determine bins
        all_vals = np.concatenate([prior_present, prior_future, post_present, post_future])
        bins = np.linspace(np.min(all_vals), np.max(all_vals), 25)

        # Plot histograms
        alpha = 0.5

        # Prior distributions
        ax.hist(prior_present, bins=bins, alpha=alpha, color=colors['prior_present'],
               label=f'Prior present', edgecolor='black', linewidth=0.3)
        ax.hist(prior_future, bins=bins, alpha=alpha, color=colors['prior_future'],
               label=f'Prior future', edgecolor='black', linewidth=0.3)

        # Posterior distributions
        ax.hist(post_present, bins=bins, alpha=alpha, color=colors['post_present'],
               label=f'Post present', edgecolor='black', linewidth=0.3)
        ax.hist(post_future, bins=bins, alpha=alpha, color=colors['post_future'],
               label=f'Post future', edgecolor='black', linewidth=0.3)

        # Add mean lines
        ymin, ymax = ax.get_ylim()
        line_height = ymax * 0.9

        if len(prior_present) > 0:
            mean_pp = np.mean(prior_present)
            ax.axvline(mean_pp, color=colors['prior_present'], linestyle='--', linewidth=2)
            ax.text(mean_pp, line_height, f'{mean_pp:.1f}', fontsize=7, color=colors['prior_present'],
                   ha='center', va='bottom')

        if len(prior_future) > 0:
            mean_pf = np.mean(prior_future)
            ax.axvline(mean_pf, color=colors['prior_future'], linestyle='--', linewidth=2)
            ax.text(mean_pf, line_height * 0.85, f'{mean_pf:.1f}', fontsize=7, color=colors['prior_future'],
                   ha='center', va='bottom')

        if len(post_present) > 0:
            mean_op = np.mean(post_present)
            ax.axvline(mean_op, color=colors['post_present'], linestyle='-', linewidth=2)
            ax.text(mean_op, line_height * 0.7, f'{mean_op:.1f}', fontsize=7, color=colors['post_present'],
                   ha='center', va='bottom')

        if len(post_future) > 0:
            mean_of = np.mean(post_future)
            ax.axvline(mean_of, color=colors['post_future'], linestyle='-', linewidth=2)
            ax.text(mean_of, line_height * 0.55, f'{mean_of:.1f}', fontsize=7, color=colors['post_future'],
                   ha='center', va='bottom')

        ax.set_xlabel('Spring Flux (m³/d)', fontsize=9)
        ax.set_ylabel('Frequency', fontsize=9)
        ax.set_title(comp['name'], fontsize=10, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=8)

        if idx == 0:
            ax.legend(fontsize=7, loc='upper right')

    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, f'{run_name}_spring_flux_distributions{file_suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved plot to: {output_path}")


if __name__ == '__main__':
    args = parse_args()

    run_name = args.run_name
    model_name = args.model_name if args.model_name else run_name

    create_spring_distribution_plot(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        spring_i=args.spring_i,
        spring_j=args.spring_j,
        suffix=args.suffix
    )
