"""
pp_budget_distributions.py - Compare prior vs posterior budget distributions

Creates distribution plots for budget differences between present and past scenarios.

Plot 1 (Differences):
- Top row: Steady state diff (kper3 - kper1)
- Middle row: Transient diff (kper4 - kper2)
- Bottom row: Seasonal change diff ((kper4-kper3) - (kper2-kper1))

Plot 2 (Raw distributions):
- Top row: Present (kper 1)
- Bottom row: Past (kper 3)

Usage:
    python scripts/pp_budget_distributions.py run_name [options]

Examples:
    python scripts/pp_budget_distributions.py run37 --run-name run37_62894 --post-iter 39
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
        description='Compare prior vs posterior budget distributions'
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


def get_iqr_limits(data, margin=0.1):
    """Get axis limits based on IQR * 1.5 or actual min/max."""
    q25 = np.percentile(data, 25)
    q75 = np.percentile(data, 75)
    iqr = q75 - q25

    iqr_min = q25 - 1.5 * iqr
    iqr_max = q75 + 1.5 * iqr

    actual_min = np.min(data)
    actual_max = np.max(data)

    # Use whichever is smaller/larger
    limit_min = max(actual_min, iqr_min)
    limit_max = min(actual_max, iqr_max)

    # Add margin
    range_val = limit_max - limit_min
    return limit_min - range_val * margin, limit_max + range_val * margin


def create_budget_distributions(run_name, model_name, post_iter=19, filter_file=None, suffix=None):
    """
    Create budget distribution plots.
    """
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"BUDGET DISTRIBUTION ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    model_ws = os.path.join('models', run_name)
    master_dir = os.path.join(model_ws, 'pest', 'master_ies')
    output_dir = os.path.join(model_ws, 'pp_budget_distributions')

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

    # Get budget types that have all needed kpers
    budget_types = ['Awanui', 'Poukawa', 'Spring', 'Confined', 'Recharge']
    valid_types = []
    for btype in budget_types:
        if btype in budget_obs:
            kpers = budget_obs[btype]
            if all(k in kpers for k in [1, 2, 3, 4]):
                valid_types.append(btype)
            else:
                print(f"  Warning: {btype} missing some kpers: {list(kpers.keys())}")

    n_cols = len(valid_types)
    if n_cols == 0:
        print("Error: No budget types have all required kpers (1, 2, 3, 4)")
        return

    print(f"\nValid budget types: {valid_types}")

    # Colors
    colors = {
        'prior': '#8d8d8d',      # grey
        'posterior': '#1f77b4',  # blue
    }

    alpha = 0.6

    # ========================================================================
    # PLOT 1: Differences
    # ========================================================================
    print("\nCreating difference plot...")

    fig, axes = plt.subplots(3, n_cols, figsize=(7.28, 7.28 * 0.6))

    # Panel labels A-O
    panel_labels = [chr(65 + i) for i in range(15)]  # A-O

    for col_idx, btype in enumerate(valid_types):
        kper_dict = budget_obs[btype]

        # Get observation names
        obs_k1 = kper_dict[1]
        obs_k2 = kper_dict[2]
        obs_k3 = kper_dict[3]
        obs_k4 = kper_dict[4]

        # --- Row 0: Steady state difference (kper 3 - kper 1) ---
        ax = axes[0, col_idx]
        panel_idx = col_idx

        if obs_k1 in post_oe.columns and obs_k3 in post_oe.columns:
            prior_diff = prior_oe[obs_k3].values - prior_oe[obs_k1].values
            post_diff = post_oe[obs_k3].values - post_oe[obs_k1].values

            xlim = get_iqr_limits(post_diff)

            bins = np.linspace(xlim[0], xlim[1], 30)
            ax.hist(prior_diff, bins=bins, alpha=alpha, color=colors['prior'],
                   label='Prior', edgecolor='black', linewidth=0.3)
            ax.hist(post_diff, bins=bins, alpha=alpha, color=colors['posterior'],
                   label='Posterior', edgecolor='black', linewidth=0.3)

            ax.set_xlim(xlim)

            ax.axvline(np.median(prior_diff), color=colors['prior'], linestyle='--', linewidth=1.5)
            ax.axvline(np.median(post_diff), color=colors['posterior'], linestyle='-', linewidth=1.5)

            ax.text(0.02, 0.98, f"Prior: {np.median(prior_diff):.1f}\nPost: {np.median(post_diff):.1f}",
                   transform=ax.transAxes, fontsize=6, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_title(f'{panel_labels[panel_idx]}: {btype}', fontsize=6, fontweight='bold', loc='left')
        ax.set_xlabel('Difference (m³/d)', fontsize=6)
        ax.set_ylabel('Frequency', fontsize=6)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3)

        if col_idx == 0:
            ax.legend(fontsize=5, loc='upper right')

        # --- Row 1: Transient difference (kper 4 - kper 2) ---
        ax = axes[1, col_idx]
        panel_idx = n_cols + col_idx

        if obs_k2 in post_oe.columns and obs_k4 in post_oe.columns:
            prior_diff = prior_oe[obs_k4].values - prior_oe[obs_k2].values
            post_diff = post_oe[obs_k4].values - post_oe[obs_k2].values

            xlim = get_iqr_limits(post_diff)

            bins = np.linspace(xlim[0], xlim[1], 30)
            ax.hist(prior_diff, bins=bins, alpha=alpha, color=colors['prior'],
                   edgecolor='black', linewidth=0.3)
            ax.hist(post_diff, bins=bins, alpha=alpha, color=colors['posterior'],
                   edgecolor='black', linewidth=0.3)

            ax.set_xlim(xlim)

            ax.axvline(np.median(prior_diff), color=colors['prior'], linestyle='--', linewidth=1.5)
            ax.axvline(np.median(post_diff), color=colors['posterior'], linestyle='-', linewidth=1.5)

            ax.text(0.02, 0.98, f"Prior: {np.median(prior_diff):.1f}\nPost: {np.median(post_diff):.1f}",
                   transform=ax.transAxes, fontsize=6, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_title(f'{panel_labels[panel_idx]}: {btype}', fontsize=6, fontweight='bold', loc='left')
        ax.set_xlabel('Difference (m³/d)', fontsize=6)
        ax.set_ylabel('Frequency', fontsize=6)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3)

        # --- Row 2: Seasonal change difference ((kper4 - kper3) - (kper2 - kper1)) ---
        ax = axes[2, col_idx]
        panel_idx = 2 * n_cols + col_idx

        if all(obs in post_oe.columns for obs in [obs_k1, obs_k2, obs_k3, obs_k4]):
            prior_past_seasonal = prior_oe[obs_k4].values - prior_oe[obs_k3].values
            prior_present_seasonal = prior_oe[obs_k2].values - prior_oe[obs_k1].values
            prior_diff = prior_past_seasonal - prior_present_seasonal

            post_past_seasonal = post_oe[obs_k4].values - post_oe[obs_k3].values
            post_present_seasonal = post_oe[obs_k2].values - post_oe[obs_k1].values
            post_diff = post_past_seasonal - post_present_seasonal

            xlim = get_iqr_limits(post_diff)

            bins = np.linspace(xlim[0], xlim[1], 30)
            ax.hist(prior_diff, bins=bins, alpha=alpha, color=colors['prior'],
                   edgecolor='black', linewidth=0.3)
            ax.hist(post_diff, bins=bins, alpha=alpha, color=colors['posterior'],
                   edgecolor='black', linewidth=0.3)

            ax.set_xlim(xlim)

            ax.axvline(np.median(prior_diff), color=colors['prior'], linestyle='--', linewidth=1.5)
            ax.axvline(np.median(post_diff), color=colors['posterior'], linestyle='-', linewidth=1.5)

            ax.text(0.02, 0.98, f"Prior: {np.median(prior_diff):.1f}\nPost: {np.median(post_diff):.1f}",
                   transform=ax.transAxes, fontsize=6, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_title(f'{panel_labels[panel_idx]}: {btype}', fontsize=6, fontweight='bold', loc='left')
        ax.set_xlabel('Difference (m³/d)', fontsize=6)
        ax.set_ylabel('Frequency', fontsize=6)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    output_path = os.path.join(output_dir, f'{run_name}_budget_differences{file_suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")

    # ========================================================================
    # PLOT 2: Raw distributions (Present kper1 and Past kper3)
    # ========================================================================
    print("\nCreating raw distribution plot...")

    fig, axes = plt.subplots(2, n_cols, figsize=(7.28, 7.28 * 0.4))

    colors_raw = {
        'prior': '#8d8d8d',
        'present': '#ff7f0e',  # orange
        'past': '#9467bd',     # purple
    }

    for col_idx, btype in enumerate(valid_types):
        kper_dict = budget_obs[btype]

        obs_k1 = kper_dict[1]
        obs_k3 = kper_dict[3]

        # --- Row 0: Present (kper 1) ---
        ax = axes[0, col_idx]

        if obs_k1 in post_oe.columns:
            prior_vals = prior_oe[obs_k1].values
            post_vals = post_oe[obs_k1].values

            xlim = get_iqr_limits(post_vals)

            bins = np.linspace(xlim[0], xlim[1], 30)
            ax.hist(prior_vals, bins=bins, alpha=alpha, color=colors_raw['prior'],
                   label='Prior', edgecolor='black', linewidth=0.3)
            ax.hist(post_vals, bins=bins, alpha=alpha, color=colors_raw['present'],
                   label='Posterior', edgecolor='black', linewidth=0.3)

            ax.set_xlim(xlim)

            ax.axvline(np.median(prior_vals), color=colors_raw['prior'], linestyle='--', linewidth=1.5)
            ax.axvline(np.median(post_vals), color=colors_raw['present'], linestyle='-', linewidth=1.5)

            ax.text(0.02, 0.98, f"Prior: {np.median(prior_vals):.1f}\nPost: {np.median(post_vals):.1f}",
                   transform=ax.transAxes, fontsize=6, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_title(f'{chr(65 + col_idx)}: {btype} (kper 1)', fontsize=6, fontweight='bold', loc='left')
        ax.set_xlabel('Budget (m³/d)', fontsize=6)
        ax.set_ylabel('Frequency', fontsize=6)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3)

        if col_idx == 0:
            ax.legend(fontsize=5, loc='upper right')

        # --- Row 1: Past (kper 3) ---
        ax = axes[1, col_idx]

        if obs_k3 in post_oe.columns:
            prior_vals = prior_oe[obs_k3].values
            post_vals = post_oe[obs_k3].values

            xlim = get_iqr_limits(post_vals)

            bins = np.linspace(xlim[0], xlim[1], 30)
            ax.hist(prior_vals, bins=bins, alpha=alpha, color=colors_raw['prior'],
                   edgecolor='black', linewidth=0.3)
            ax.hist(post_vals, bins=bins, alpha=alpha, color=colors_raw['past'],
                   edgecolor='black', linewidth=0.3)

            ax.set_xlim(xlim)

            ax.axvline(np.median(prior_vals), color=colors_raw['prior'], linestyle='--', linewidth=1.5)
            ax.axvline(np.median(post_vals), color=colors_raw['past'], linestyle='-', linewidth=1.5)

            ax.text(0.02, 0.98, f"Prior: {np.median(prior_vals):.1f}\nPost: {np.median(post_vals):.1f}",
                   transform=ax.transAxes, fontsize=6, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

        ax.set_title(f'{chr(65 + n_cols + col_idx)}: {btype} (kper 3)', fontsize=6, fontweight='bold', loc='left')
        ax.set_xlabel('Budget (m³/d)', fontsize=6)
        ax.set_ylabel('Frequency', fontsize=6)
        ax.tick_params(labelsize=6)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()

    output_path = os.path.join(output_dir, f'{run_name}_budget_raw{file_suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")

    # Save summary statistics
    print("\nSaving summary statistics...")
    summary_rows = []

    for btype in valid_types:
        kper_dict = budget_obs[btype]

        for kper in [1, 2, 3, 4]:
            obs_name = kper_dict[kper]
            if obs_name in post_oe.columns:
                prior_vals = prior_oe[obs_name].values
                post_vals = post_oe[obs_name].values

                summary_rows.append({
                    'budget_type': btype,
                    'kper': kper,
                    'prior_median': np.median(prior_vals),
                    'prior_mean': np.mean(prior_vals),
                    'prior_std': np.std(prior_vals),
                    'post_median': np.median(post_vals),
                    'post_mean': np.mean(post_vals),
                    'post_std': np.std(post_vals),
                })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, f'{run_name}_budget_summary{file_suffix}.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"  Saved: {summary_path}")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_budget_distributions(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        suffix=args.suffix
    )
