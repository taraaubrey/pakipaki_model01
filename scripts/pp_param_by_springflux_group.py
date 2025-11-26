"""
pp_param_by_springflux_group.py - Analyze parameter distributions by spring flux change groups

Groups realizations based on whether steady state spring flux increased or decreased
from present to past, then compares parameter distributions between groups.

Usage:
    python scripts/pp_param_by_springflux_group.py run_name [options]

Examples:
    python scripts/pp_param_by_springflux_group.py run37 --run-name run37_62894 --post-iter 39 --suffix _fh18
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
        description='Analyze parameter distributions by spring flux change groups'
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
    parser.add_argument('--post-iter', type=int, default=39,
                       help='Posterior iteration (default: 19)')
    parser.add_argument('--suffix', type=str, default='',
                       help='Suffix for input data file (e.g., "_fh18")')
    return parser.parse_args()


def analyze_param_by_groups(run_name, model_name, post_iter=19, suffix=''):
    """
    Analyze parameter distributions by spring flux change groups.
    """
    print(f"\n{'='*80}")
    print(f"PARAMETER ANALYSIS BY SPRING FLUX GROUPS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    data_dir = os.path.join('models', run_name, 'pp_predictions_springflux', 'data')
    output_dir = os.path.join('models', run_name, 'pp_param_by_springflux_group')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load spring flux data
    flux_file = os.path.join(data_dir, f'{run_name}_spring_flux_data{suffix}.csv')
    print(f"\nLoading spring flux data from: {flux_file}")
    flux_df = pd.read_csv(flux_file)

    # Filter for steady state, posterior
    steady_post = flux_df[(flux_df['comparison'] == 'Steady state') &
                          (flux_df['iteration'] == 'posterior')]

    # Pivot to get present and past values per realization
    present = steady_post[steady_post['period'] == 'present'].set_index('realization')['value']
    past = steady_post[steady_post['period'] == 'past'].set_index('realization')['value']

    # Calculate difference (past - present)
    # Positive = past flux more negative (more outflow in past)
    # Negative = present flux more negative (more outflow in present)
    diff = past - present

    # Group realizations
    # Note: spring flux is negative (outflow), so:
    # - If past is more negative than present, diff < 0 (past had more outflow)
    # - If present is more negative than past, diff > 0 (present has more outflow)
    group_pos = diff[diff > 0].index.tolist()  # Present has more outflow
    group_neg = diff[diff < 0].index.tolist()  # Past had more outflow

    print(f"\nSteady state spring flux analysis:")
    print(f"  Total realizations: {len(diff)}")
    print(f"  Group 1 (diff > 0, present more outflow): {len(group_pos)} realizations")
    print(f"  Group 2 (diff < 0, past more outflow): {len(group_neg)} realizations")

    # Print some statistics
    print(f"\nDifference statistics (past - present):")
    print(f"  Min: {diff.min():.2f}")
    print(f"  Max: {diff.max():.2f}")
    print(f"  Mean: {diff.mean():.2f}")
    print(f"  Median: {diff.median():.2f}")

    # Save group assignments
    group_df = pd.DataFrame({
        'realization': diff.index,
        'present_flux': present.values,
        'past_flux': past.values,
        'diff': diff.values,
        'group': ['present_more' if d > 0 else 'past_more' for d in diff.values]
    })
    group_file = os.path.join(output_dir, f'{run_name}_springflux_groups{suffix}.csv')
    group_df.to_csv(group_file, index=False)
    print(f"\nSaved group assignments to: {group_file}")

    # Load PST to get parameter info
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    print(f"\nLoading PST: {pst_file}")
    pst = pyemu.Pst(pst_file)
    par_data = pst.parameter_data

    # Get parameter groups
    par_groups = par_data['pargp'].unique()
    print(f"Parameter groups: {len(par_groups)}")

    # Load posterior parameter ensemble
    par_file = os.path.join(master_dir, f'{model_name}.{post_iter}.par.csv')
    print(f"Loading parameter ensemble: {par_file}")
    par_oe = pd.read_csv(par_file, index_col=0)

    # Filter to only realizations in our groups
    # Handle mixed types by converting to strings for comparison
    all_reals = group_pos + group_neg
    par_index_str = [str(i) for i in par_oe.index]
    all_reals_str = [str(r) for r in all_reals]

    mask = [i for i, idx in enumerate(par_index_str) if idx in all_reals_str]
    par_oe = par_oe.iloc[mask]

    # Also convert group lists to match par_oe index type
    par_oe_index_set = set(str(i) for i in par_oe.index)
    group_pos = [r for r in group_pos if str(r) in par_oe_index_set]
    group_neg = [r for r in group_neg if str(r) in par_oe_index_set]

    print(f"  Loaded {len(par_oe)} realizations")

    # Get parameter group statistics for each group
    print("\nCalculating parameter statistics by group...")

    stats_rows = []

    for pargp in sorted(par_groups):
        # Get parameters in this group
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
        pars_in_group = [p for p in pars_in_group if p in par_oe.columns]

        if len(pars_in_group) == 0:
            continue

        # Check if log-transformed
        partrans = par_data.loc[pars_in_group[0], 'partrans'] if len(pars_in_group) > 0 else 'none'
        is_log = partrans == 'log'

        # Get values for each group - use string matching for index
        pos_reals = [r for r in group_pos if str(r) in par_oe_index_set]
        neg_reals = [r for r in group_neg if str(r) in par_oe_index_set]

        if len(pos_reals) == 0 or len(neg_reals) == 0:
            continue

        # Calculate group statistics (mean across all parameters in group)
        # Convert realization IDs to match index type
        pos_idx = [par_oe.index[par_index_str.index(str(r))] for r in pos_reals if str(r) in par_index_str]
        neg_idx = [par_oe.index[par_index_str.index(str(r))] for r in neg_reals if str(r) in par_index_str]

        pos_vals = par_oe.loc[pos_idx, pars_in_group].values.flatten()
        neg_vals = par_oe.loc[neg_idx, pars_in_group].values.flatten()

        # If log-transformed, convert to real space for statistics
        if is_log:
            pos_vals_real = 10 ** pos_vals
            neg_vals_real = 10 ** neg_vals
        else:
            pos_vals_real = pos_vals
            neg_vals_real = neg_vals

        stats_rows.append({
            'pargp': pargp,
            'n_pars': len(pars_in_group),
            'partrans': partrans,
            'pos_n': len(pos_reals),
            'pos_median': np.median(pos_vals_real),
            'pos_mean': np.mean(pos_vals_real),
            'pos_std': np.std(pos_vals_real),
            'pos_p10': np.percentile(pos_vals_real, 10),
            'pos_p90': np.percentile(pos_vals_real, 90),
            'neg_n': len(neg_reals),
            'neg_median': np.median(neg_vals_real),
            'neg_mean': np.mean(neg_vals_real),
            'neg_std': np.std(neg_vals_real),
            'neg_p10': np.percentile(neg_vals_real, 10),
            'neg_p90': np.percentile(neg_vals_real, 90),
        })

    stats_df = pd.DataFrame(stats_rows)
    stats_file = os.path.join(output_dir, f'{run_name}_param_stats_by_group{suffix}.csv')
    stats_df.to_csv(stats_file, index=False)
    print(f"Saved parameter statistics to: {stats_file}")

    # Print summary table
    print("\n" + "="*100)
    print("PARAMETER GROUP STATISTICS BY SPRING FLUX GROUP")
    print("="*100)
    print(f"{'Parameter Group':<25} {'Trans':<6} {'Present>Past Median':<20} {'Past>Present Median':<20}")
    print("-"*100)
    for _, row in stats_df.iterrows():
        print(f"{row['pargp']:<25} {row['partrans']:<6} {row['pos_median']:<20.4g} {row['neg_median']:<20.4g}")

    # Create distribution plots
    print("\nCreating distribution plots...")

    # Determine number of parameter groups to plot
    n_groups = len(stats_df)
    n_cols = 4
    n_rows = int(np.ceil(n_groups / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(7.28, n_rows * 1.5))
    axes = axes.flatten() if n_groups > 1 else [axes]

    colors = {
        'pos': '#ff7f0e',  # orange - present more outflow
        'neg': '#9467bd',  # purple - past more outflow
    }

    for idx, (_, row) in enumerate(stats_df.iterrows()):
        if idx >= len(axes):
            break

        ax = axes[idx]
        pargp = row['pargp']

        # Get parameters in this group
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
        pars_in_group = [p for p in pars_in_group if p in par_oe.columns]

        if len(pars_in_group) == 0:
            ax.set_visible(False)
            continue

        # Check if log-transformed
        partrans = par_data.loc[pars_in_group[0], 'partrans']
        is_log = partrans == 'log'

        # Get values for each group - use string matching for index
        pos_idx = [par_oe.index[par_index_str.index(str(r))] for r in group_pos if str(r) in par_index_str]
        neg_idx = [par_oe.index[par_index_str.index(str(r))] for r in group_neg if str(r) in par_index_str]

        if len(pos_idx) == 0 or len(neg_idx) == 0:
            ax.set_visible(False)
            continue

        pos_vals = par_oe.loc[pos_idx, pars_in_group].values.flatten()
        neg_vals = par_oe.loc[neg_idx, pars_in_group].values.flatten()

        # Convert to real space if log-transformed
        if is_log:
            pos_vals = 10 ** pos_vals
            neg_vals = 10 ** neg_vals

        # Filter out NaN and Inf values
        pos_vals = pos_vals[np.isfinite(pos_vals)]
        neg_vals = neg_vals[np.isfinite(neg_vals)]

        if len(pos_vals) == 0 or len(neg_vals) == 0:
            ax.text(0.5, 0.5, 'Invalid values', ha='center', va='center', fontsize=6)
            ax.set_title(f"{chr(65+idx)}: {pargp}", fontsize=6, fontweight='bold', loc='left')
            continue

        # Plot histograms
        # Use common bins based on combined data
        all_vals = np.concatenate([pos_vals, neg_vals])

        # Use IQR-based limits
        q25, q75 = np.percentile(all_vals, [25, 75])
        iqr = q75 - q25
        x_min = max(np.min(all_vals), q25 - 1.5 * iqr)
        x_max = min(np.max(all_vals), q75 + 1.5 * iqr)

        # Safeguard against NaN/Inf
        if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min >= x_max:
            ax.text(0.5, 0.5, 'Invalid range', ha='center', va='center', fontsize=6)
            ax.set_title(f"{chr(65+idx)}: {pargp}", fontsize=6, fontweight='bold', loc='left')
            continue

        bins = np.linspace(x_min, x_max, 30)

        ax.hist(pos_vals, bins=bins, alpha=0.6, color=colors['pos'],
               label='Present>Past', edgecolor='black', linewidth=0.3, density=True)
        ax.hist(neg_vals, bins=bins, alpha=0.6, color=colors['neg'],
               label='Past>Present', edgecolor='black', linewidth=0.3, density=True)

        # Add median lines
        ax.axvline(np.median(pos_vals), color=colors['pos'], linestyle='-', linewidth=1.5)
        ax.axvline(np.median(neg_vals), color=colors['neg'], linestyle='-', linewidth=1.5)

        ax.set_xlim(x_min, x_max)
        ax.set_title(f"{chr(65+idx)}: {pargp}", fontsize=6, fontweight='bold', loc='left')
        ax.tick_params(labelsize=5)
        ax.grid(True, alpha=0.3)

        if idx == 0:
            ax.legend(fontsize=5, loc='upper right')

    # Hide unused axes
    for idx in range(len(stats_df), len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()

    # Save plot
    plot_file = os.path.join(output_dir, f'{run_name}_param_distributions_by_group{suffix}.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved distribution plot to: {plot_file}")

    # Create a more detailed plot for key parameter groups
    key_groups = ['k', 'ss', 'rch', 'ghbcond', 'ghbhead']
    key_pargps = []
    for kg in key_groups:
        for pargp in par_groups:
            if kg in pargp.lower():
                key_pargps.append(pargp)

    if len(key_pargps) > 0:
        print(f"\nCreating detailed plots for {len(key_pargps)} key parameter groups...")

        n_key = len(key_pargps)
        n_cols_key = min(4, n_key)
        n_rows_key = int(np.ceil(n_key / n_cols_key))

        fig, axes = plt.subplots(n_rows_key, n_cols_key, figsize=(7.28, n_rows_key * 2))
        if n_key == 1:
            axes = [axes]
        else:
            axes = axes.flatten()

        for idx, pargp in enumerate(key_pargps):
            if idx >= len(axes):
                break

            ax = axes[idx]

            # Get parameters in this group
            pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
            pars_in_group = [p for p in pars_in_group if p in par_oe.columns]

            if len(pars_in_group) == 0:
                ax.set_visible(False)
                continue

            # Check if log-transformed
            partrans = par_data.loc[pars_in_group[0], 'partrans']
            is_log = partrans == 'log'

            # Get values for each group - use string matching for index
            pos_idx = [par_oe.index[par_index_str.index(str(r))] for r in group_pos if str(r) in par_index_str]
            neg_idx = [par_oe.index[par_index_str.index(str(r))] for r in group_neg if str(r) in par_index_str]

            if len(pos_idx) == 0 or len(neg_idx) == 0:
                ax.set_visible(False)
                continue

            pos_vals = par_oe.loc[pos_idx, pars_in_group].values.flatten()
            neg_vals = par_oe.loc[neg_idx, pars_in_group].values.flatten()

            # Convert to real space if log-transformed
            if is_log:
                pos_vals = 10 ** pos_vals
                neg_vals = 10 ** neg_vals

            # Filter out NaN and Inf values
            pos_vals = pos_vals[np.isfinite(pos_vals)]
            neg_vals = neg_vals[np.isfinite(neg_vals)]

            if len(pos_vals) == 0 or len(neg_vals) == 0:
                ax.text(0.5, 0.5, 'Invalid values', ha='center', va='center', fontsize=6)
                ax.set_title(f"{chr(65+idx)}: {pargp}", fontsize=6, fontweight='bold', loc='left')
                continue

            # Box plots
            data = [pos_vals, neg_vals]
            bp = ax.boxplot(data, labels=['Present>Past', 'Past>Present'],
                           patch_artist=True, showfliers=False)

            bp['boxes'][0].set_facecolor(colors['pos'])
            bp['boxes'][1].set_facecolor(colors['neg'])

            for box in bp['boxes']:
                box.set_alpha(0.6)

            ax.set_title(f"{chr(65+idx)}: {pargp}", fontsize=6, fontweight='bold', loc='left')
            ax.tick_params(labelsize=5)
            ax.grid(True, alpha=0.3, axis='y')

            # Add median text
            pos_med = np.median(pos_vals)
            neg_med = np.median(neg_vals)
            ax.text(0.02, 0.98, f"Med: {pos_med:.3g} | {neg_med:.3g}",
                   transform=ax.transAxes, fontsize=5, verticalalignment='top')

        # Hide unused axes
        for idx in range(len(key_pargps), len(axes)):
            axes[idx].set_visible(False)

        plt.tight_layout()

        # Save plot
        box_file = os.path.join(output_dir, f'{run_name}_param_boxplots_by_group{suffix}.png')
        plt.savefig(box_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved boxplot to: {box_file}")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    analyze_param_by_groups(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        suffix=args.suffix
    )
