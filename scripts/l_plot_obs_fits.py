"""
Plot observation fits for filtered posterior ensemble.

This script creates comprehensive observation fit plots including:
- 1:1 plots with boxplots
- CV reduction plots
- Distribution histograms and time series

Usage:
    python scripts/l_plot_obs_fits.py local_run32 --iteration 16 --arr-h-threshold 18
    python scripts/l_plot_obs_fits.py local_run32 -i 16 --arr-h-threshold 300
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyemu

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from setup import MODEL_NAME as DEFAULT_MODEL_NAME


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Plot observation fits for filtered posterior ensemble',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default=DEFAULT_MODEL_NAME,
        help=f'Model/run name (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )

    parser.add_argument(
        '--iteration', '-i',
        type=int,
        required=True,
        help='Posterior iteration to analyze'
    )

    parser.add_argument(
        '--arr-h-threshold',
        type=float,
        default=18.0,
        help='Minimum mean arr-h threshold for filtering realizations (default: 18)'
    )

    return parser.parse_args()


def filter_ensemble(oe, obs_data, arr_h_threshold):
    """
    Filter ensemble based on criteria:
    1. Mean arr-h > threshold
    2. Confined q budget not negative

    Returns filtered observation ensemble and list of kept realization IDs.
    """
    print(f"\nFiltering ensemble...")
    print(f"  Initial ensemble size: {oe.shape[0]} realizations")

    # Get arr-h observation names
    arr_h_obs = obs_data[obs_data['oname'] == 'arr-h'].index.tolist()

    # Get confined budget observation names
    conf_budget_obs = obs_data[obs_data['oname'].str.contains('budget.*confined', case=False, na=False)].index.tolist()

    # Filter 1: Mean arr-h > threshold
    if len(arr_h_obs) > 0:
        arr_h_cols = [col for col in oe.columns if col in arr_h_obs]
        if len(arr_h_cols) > 0:
            mean_arr_h = oe[arr_h_cols].mean(axis=1)
            mask_arr_h = mean_arr_h > arr_h_threshold
            print(f"  Filter 1 (mean arr-h > {arr_h_threshold}): {mask_arr_h.sum()} realizations pass")
        else:
            mask_arr_h = pd.Series([True] * len(oe), index=oe.index)
            print(f"  Filter 1 (mean arr-h): No arr-h observations found, skipping")
    else:
        mask_arr_h = pd.Series([True] * len(oe), index=oe.index)
        print(f"  Filter 1 (mean arr-h): No arr-h observations found, skipping")

    # Filter 2: Confined budget not negative
    if len(conf_budget_obs) > 0:
        conf_budget_cols = [col for col in oe.columns if col in conf_budget_obs]
        if len(conf_budget_cols) > 0:
            # Check if ANY confined budget obs is negative
            mask_budget = (oe[conf_budget_cols] >= 0).all(axis=1)
            print(f"  Filter 2 (confined budget >= 0): {mask_budget.sum()} realizations pass")
        else:
            mask_budget = pd.Series([True] * len(oe), index=oe.index)
            print(f"  Filter 2 (confined budget): No confined budget observations found, skipping")
    else:
        mask_budget = pd.Series([True] * len(oe), index=oe.index)
        print(f"  Filter 2 (confined budget): No confined budget observations found, skipping")

    # Combine filters
    mask_combined = mask_arr_h & mask_budget

    oe_filtered = oe[mask_combined]
    kept_ids = oe_filtered.index.tolist()

    print(f"  Final filtered ensemble size: {len(kept_ids)} realizations")
    print(f"  Dropped: {len(oe) - len(kept_ids)} realizations")

    return oe_filtered, kept_ids


def get_non_inequality_obs(obs_data, phi_factors):
    """
    Get list of non-zero weighted observations that are not inequality constraints.
    Excludes observations with groups starting with 'less_' or 'greater_'.
    """
    # Get non-zero weighted obs
    non_zero = phi_factors[phi_factors['weight'] > 0]['group'].tolist()

    # Filter out inequality constraints
    non_inequality = [obs for obs in non_zero
                      if not obs.startswith('less_') and not obs.startswith('greater_')]

    # Get observation names from obs_data
    obs_list = []
    for group in non_inequality:
        if group in obs_data.index:
            obs_list.append(group)

    return obs_list


def load_obs_noise(model_ws):
    """Load observation + noise data."""
    obs_noise_file = os.path.join(model_ws, 'truth', 'output.obs+noise.csv')
    if os.path.exists(obs_noise_file):
        return pd.read_csv(obs_noise_file, index_col=0)
    else:
        print(f"Warning: obs+noise file not found: {obs_noise_file}")
        return None


def load_truth(model_ws):
    """Load all truth observation data and combine into single DataFrame."""
    # Try different truth file names
    truth_files = [
        'output.sample_heads.truth.csv',
        'output.budget.truth.csv',
        'output.sample_recession_rates.truth.csv',
        'output.GHB_AW_fluxes.truth.csv'
    ]

    all_truth = []
    for fname in truth_files:
        truth_file = os.path.join(model_ws, 'truth', fname)
        if os.path.exists(truth_file):
            try:
                df = pd.read_csv(truth_file, index_col=0)
                all_truth.append(df)
                print(f"  Loaded truth: {fname} ({len(df)} rows)")
            except Exception as e:
                print(f"  Warning: Could not load {fname}: {e}")

    if len(all_truth) == 0:
        print(f"Warning: No truth files found in {os.path.join(model_ws, 'truth')}")
        return None

    # Combine all truth data
    combined = pd.concat(all_truth)
    # Handle duplicate indices by keeping first
    combined = combined[~combined.index.duplicated(keep='first')]
    print(f"  Combined truth data: {len(combined)} observations")
    return combined


def main():
    # args = parse_args()
    # model_name = args.model_name
    # iteration = args.iteration
    # arr_h_threshold = args.arr_h_threshold
    model_name = 'local_run34'
    iteration = 7
    arr_h_threshold = 18.0

    print(f"\n{'='*80}")
    print(f"OBSERVATION FIT ANALYSIS: {model_name}")
    print(f"{'='*80}")
    print(f"Posterior iteration: {iteration}")
    print(f"arr-h threshold: {arr_h_threshold}")

    # Setup paths
    model_ws = os.path.join('models', model_name)
    master_dir = os.path.join(model_ws, 'pest', 'master_ies')
    output_dir = os.path.join(model_ws, 'figures')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST and observation data
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Load phi factors
    phi_factors_file = os.path.join(master_dir, 'phi_factors.csv')
    phi_factors = pd.read_csv(phi_factors_file, header=None, names=['group', 'weight'])

    # Load observation ensembles
    print(f"\nLoading observation ensembles...")
    prior_oe_file = os.path.join(master_dir, f'{model_name}.0.obs.csv')
    post_oe_file = os.path.join(master_dir, f'{model_name}.{iteration}.obs.csv')

    prior_oe = pd.read_csv(prior_oe_file, index_col=0)
    post_oe = pd.read_csv(post_oe_file, index_col=0)

    print(f"  Prior ensemble: {prior_oe.shape}")
    print(f"  Posterior ensemble: {post_oe.shape}")

    # Filter posterior ensemble
    post_oe_filtered, kept_ids = filter_ensemble(post_oe, obs_data, arr_h_threshold)

    # Load truth and obs+noise
    # truth = load_truth(model_ws)
    obs_noise_file = os.path.join(master_dir, f'{model_name}.obs+noise.csv')
    obs_noise = pd.read_csv(obs_noise_file, index_col=0)

    # Get non-inequality observations
    non_ineq_obs = get_non_inequality_obs(obs_data, phi_factors)
    print(f"\nNon-inequality observations with non-zero weights: {len(non_ineq_obs)}")

    # Create figure
    print(f"\nCreating plots...")
    width = 12
    height = 0.8 * width
    fig = plt.figure(figsize=(width, height))
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.35)

    # Plot 1: 1:1 plot with boxplots [0:2, 0:2]
    ax1 = fig.add_subplot(gs[0:2, 0:2])
    plot_1to1_boxplot(ax1, prior_oe, post_oe_filtered, obs_data, non_ineq_obs, truth, obs_noise)

    # Plot 2: CV reduction [0:2, 2]
    ax2 = fig.add_subplot(gs[0:2, 2])
    plot_cv_reduction(ax2, prior_oe, post_oe_filtered, non_ineq_obs)

    # Plot 3: kper 1 histogram [2, 0]
    ax3 = fig.add_subplot(gs[2, 0])
    plot_kper1_histogram(ax3, prior_oe, post_oe_filtered, obs_data, truth, obs_noise, kper=1)

    # Plot 4: kper 2 time series [2, 1:3]
    ax4 = fig.add_subplot(gs[2, 1:3])
    plot_kper2_timeseries(ax4, prior_oe, post_oe_filtered, obs_data, truth, obs_noise, kper=2)

    # Save figure
    fig_path = os.path.join(output_dir, f'{model_name}_obs_fits_iter{iteration}.png')
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"\nFigure saved to: {fig_path}")

    print(f"\n{'='*80}")
    print("ANALYSIS COMPLETE!")
    print(f"{'='*80}\n")


def plot_1to1_boxplot(ax, prior_oe, post_oe, obs_data, obs_list, truth, obs_noise):
    """Plot 1:1 comparison with boxplots for ensemble."""
    ax.set_title('A. Observed vs simulated (1:1)', fontsize=8, loc='left', fontweight='bold')

    # Get observed values from obs_data
    obs_values = []
    sim_data = []

    for obs in obs_list:
        if obs in obs_data.index and obs in post_oe.columns:
            obs_val = obs_data.loc[obs, 'obsval']
            obs_values.append(obs_val)
            sim_data.append(post_oe[obs].values)

    if len(obs_values) == 0:
        ax.text(0.5, 0.5, 'No observations available', ha='center', va='center',
                transform=ax.transAxes, fontsize=8)
        return

    # Create boxplot for each observation
    positions = np.array(obs_values)
    bp = ax.boxplot(sim_data, positions=positions, widths=np.ptp(positions)/len(positions)/3,
                     vert=False, patch_artist=True, manage_ticks=False,
                     boxprops=dict(facecolor='lightblue', alpha=0.6),
                     medianprops=dict(color='darkblue', linewidth=1.5),
                     whiskerprops=dict(color='blue'),
                     capprops=dict(color='blue'))

    # Add 1:1 line
    lims = [min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1])]
    ax.plot(lims, lims, 'k--', alpha=0.5, zorder=0, label='1:1 line')

    ax.set_xlabel('Observed', fontsize=8)
    ax.set_ylabel('Simulated', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3, zorder=0)
    ax.set_axisbelow(True)
    ax.legend(fontsize=6, loc='best')


def plot_cv_reduction(ax, prior_oe, post_oe, obs_list):
    """Plot CV reduction from prior to posterior."""
    ax.set_title('B. CV reduction', fontsize=8, loc='left', fontweight='bold')

    # Calculate CV for prior and posterior
    # CV = std / mean * 100
    cv_reductions = []
    obs_labels = []

    for obs in obs_list:
        if obs in prior_oe.columns and obs in post_oe.columns:
            # Prior CV
            prior_mean = prior_oe[obs].mean()
            prior_std = prior_oe[obs].std()
            cv_prior = (prior_std / abs(prior_mean) * 100) if abs(prior_mean) > 1e-10 else 0

            # Posterior CV
            post_mean = post_oe[obs].mean()
            post_std = post_oe[obs].std()
            cv_post = (post_std / abs(post_mean) * 100) if abs(post_mean) > 1e-10 else 0

            # CV reduction (percentage)
            cv_reduction = ((cv_prior - cv_post) / cv_prior * 100) if cv_prior > 0 else 0

            cv_reductions.append(cv_reduction)
            obs_labels.append(obs)

    if len(cv_reductions) == 0:
        ax.text(0.5, 0.5, 'No observations available', ha='center', va='center',
                transform=ax.transAxes, fontsize=8)
        return

    # Sort by CV reduction
    sorted_indices = np.argsort(cv_reductions)
    cv_reductions = np.array(cv_reductions)[sorted_indices]
    obs_labels = np.array(obs_labels)[sorted_indices]

    # Plot horizontal bar chart
    y_pos = np.arange(len(cv_reductions))
    colors = ['green' if cv > 0 else 'red' for cv in cv_reductions]

    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.3, axis='x', zorder=0)
    ax.barh(y_pos, cv_reductions, color=colors, alpha=0.7, edgecolor='black', zorder=3)

    ax.set_yticks(y_pos[::max(1, len(y_pos)//10)])  # Show max 10 labels
    ax.set_yticklabels(obs_labels[::max(1, len(y_pos)//10)], fontsize=6)
    ax.set_xlabel('CV reduction (%)', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.axvline(x=0, color='black', linewidth=1, linestyle='--', alpha=0.5)


def plot_kper1_histogram(ax, prior_oe, post_oe, obs_data, truth, obs_noise, kper=1):
    """Plot histogram for kper 1 observations."""
    ax.set_title(f'C. Kper {kper} distribution', fontsize=8, loc='left', fontweight='bold')

    # Get observations for this kper
    kper_obs = obs_data[obs_data['kper'] == kper].index.tolist()
    kper_obs = [obs for obs in kper_obs if obs in post_oe.columns]

    if len(kper_obs) == 0:
        ax.text(0.5, 0.5, f'No observations for kper {kper}', ha='center', va='center',
                transform=ax.transAxes, fontsize=8)
        return

    # Collect all values from all observations for this kper
    prior_vals = []
    post_vals = []
    truth_vals = []
    obs_noise_vals = []

    for obs in kper_obs:
        if obs in prior_oe.columns:
            prior_vals.extend(prior_oe[obs].values)
        if obs in post_oe.columns:
            post_vals.extend(post_oe[obs].values)
        if truth is not None and obs in truth.index:
            truth_vals.append(truth.loc[obs])
        if obs_noise is not None and obs in obs_noise.index:
            obs_noise_vals.append(obs_noise.loc[obs])

    # Plot histograms
    if len(prior_vals) > 0:
        ax.hist(prior_vals, bins=30, alpha=0.5, color='blue', label='Prior', density=True)
    if len(post_vals) > 0:
        ax.hist(post_vals, bins=30, alpha=0.5, color='red', label='Posterior', density=True)

    # Plot truth and obs+noise as vertical lines
    if len(truth_vals) > 0:
        for val in truth_vals:
            ax.axvline(val, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.plot([], [], 'k--', linewidth=1.5, label='Truth')

    if len(obs_noise_vals) > 0:
        for val in obs_noise_vals:
            ax.axvline(val, color='green', linestyle=':', linewidth=1.5, alpha=0.7)
        ax.plot([], [], 'g:', linewidth=1.5, label='Obs+noise')

    ax.set_xlabel('Value', fontsize=8)
    ax.set_ylabel('Density', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=6, loc='best')


def plot_kper2_timeseries(ax, prior_oe, post_oe, obs_data, truth, obs_noise, kper=2):
    """Plot time series for kper 2 observations with 5-95th percentile envelopes."""
    ax.set_title(f'D. Kper {kper} time series', fontsize=8, loc='left', fontweight='bold')

    # Get observations for this kper, grouped by kstp
    kper_obs = obs_data[obs_data['kper'] == kper]

    if len(kper_obs) == 0:
        ax.text(0.5, 0.5, f'No observations for kper {kper}', ha='center', va='center',
                transform=ax.transAxes, fontsize=8)
        return

    # Get unique time steps
    kstps = sorted(kper_obs['kstp'].unique())

    # For each kstp, calculate percentiles
    prior_p05 = []
    prior_p50 = []
    prior_p95 = []
    post_p05 = []
    post_p50 = []
    post_p95 = []
    truth_vals = []
    obs_noise_p05 = []
    obs_noise_p95 = []

    for kstp in kstps:
        kstp_obs = kper_obs[kper_obs['kstp'] == kstp].index.tolist()
        kstp_obs = [obs for obs in kstp_obs if obs in post_oe.columns]

        if len(kstp_obs) > 0:
            # Collect all values for this time step
            prior_all = []
            post_all = []

            for obs in kstp_obs:
                if obs in prior_oe.columns:
                    prior_all.extend(prior_oe[obs].values)
                if obs in post_oe.columns:
                    post_all.extend(post_oe[obs].values)

            # Calculate percentiles
            if len(prior_all) > 0:
                prior_p05.append(np.percentile(prior_all, 5))
                prior_p50.append(np.percentile(prior_all, 50))
                prior_p95.append(np.percentile(prior_all, 95))
            else:
                prior_p05.append(np.nan)
                prior_p50.append(np.nan)
                prior_p95.append(np.nan)

            if len(post_all) > 0:
                post_p05.append(np.percentile(post_all, 5))
                post_p50.append(np.percentile(post_all, 50))
                post_p95.append(np.percentile(post_all, 95))
            else:
                post_p05.append(np.nan)
                post_p50.append(np.nan)
                post_p95.append(np.nan)

            # Get truth and obs+noise for this time step (take first obs)
            if truth is not None and kstp_obs[0] in truth.index:
                truth_vals.append(truth.loc[kstp_obs[0]])
            else:
                truth_vals.append(np.nan)

    # Plot envelopes
    if len(kstps) > 0:
        ax.fill_between(kstps, prior_p05, prior_p95, alpha=0.2, color='blue', label='Prior 5-95%')
        ax.plot(kstps, prior_p50, 'b-', linewidth=1.5, label='Prior median')

        ax.fill_between(kstps, post_p05, post_p95, alpha=0.2, color='red', label='Posterior 5-95%')
        ax.plot(kstps, post_p50, 'r-', linewidth=1.5, label='Posterior median')

        # Plot truth
        if len(truth_vals) > 0 and not all(np.isnan(truth_vals)):
            ax.plot(kstps, truth_vals, 'k--', linewidth=2, label='Truth')

    ax.set_xlabel('Time step', fontsize=8)
    ax.set_ylabel('Value', fontsize=8)
    ax.tick_params(labelsize=8)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=6, loc='best')


if __name__ == '__main__':
    sys.exit(main())
