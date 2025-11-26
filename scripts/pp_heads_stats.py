"""
pp_heads_stats.py - Extract head array statistics from posterior ensemble

Extracts arr-h observations for a specific kper/kstp and calculates statistics
(mean, median, percentiles) across realizations for each cell.

Usage:
    python scripts/pp_heads_stats.py run_name [options]

Examples:
    python scripts/pp_heads_stats.py local_run34 --kper 1 --kstp 1 --post-iter 19
    python scripts/pp_heads_stats.py local_run34 --kper 2 --kstp 52 --post-iter 19 --filter-file models/local_run34/filtered_realizations.csv
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
        description='Extract head array statistics from posterior ensemble'
    )
    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default='run36',
        help=f'Model name for file prefixes (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )
    parser.add_argument('--run-name', '-r', type=str, default='run36_62859_noGHBheads',
                       help='Run name for directory path (default: same as model_name)')
    parser.add_argument('--post-iter', type=int, default=39,
                       help='Posterior iteration (default: 39)')
    parser.add_argument('--kper', type=int, default=1,
                       help='Stress period (default: 1)')
    parser.add_argument('--kstp', type=int, default=1,
                       help='Time step (default: 1)')
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


def create_heads_stats(run_name, model_name, post_iter=19, kper=1, kstp=1,
                       filter_file=None, suffix=None):
    """
    Extract head array statistics from posterior ensemble.
    """
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"HEAD ARRAY STATISTICS: {run_name}")
    print(f"kper={kper}, kstp={kstp}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_heads')
    arrays_dir = os.path.join(output_dir, 'arrays')

    for d in [output_dir, arrays_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data.copy()

    # Filter to arr-h observations for this kper/kstp
    arr_h = obs_data[obs_data['oname'] == 'arr-h'].copy()
    print(f"\nTotal arr-h observations: {len(arr_h)}")

    # Convert kper/kstp to numeric
    arr_h['kper'] = pd.to_numeric(arr_h['kper'], errors='coerce')
    arr_h['kstp'] = pd.to_numeric(arr_h['kstp'], errors='coerce')

    # Filter by kper/kstp
    arr_h_filtered = arr_h[(arr_h['kper'] == kper) & (arr_h['kstp'] == kstp)]
    print(f"arr-h observations for kper={kper}, kstp={kstp}: {len(arr_h_filtered)}")

    if len(arr_h_filtered) == 0:
        print("Error: No observations found for this kper/kstp")
        return

    # Get grid dimensions from i, j indices
    arr_h_filtered['i'] = pd.to_numeric(arr_h_filtered['i'], errors='coerce')
    arr_h_filtered['j'] = pd.to_numeric(arr_h_filtered['j'], errors='coerce')

    nrow = int(arr_h_filtered['i'].max()) + 1
    ncol = int(arr_h_filtered['j'].max()) + 1
    print(f"Grid dimensions: {nrow} rows x {ncol} cols")

    # Load filtered realizations
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations")

    # Load posterior observation ensemble
    print(f"\nLoading posterior observation ensemble (iter {post_iter})...")
    post_file = os.path.join(master_dir, f'{model_name}.{post_iter}.obs.csv')

    # Get observation names for this kper/kstp
    obs_names = arr_h_filtered.index.tolist()

    # Read only needed columns
    post_header = pd.read_csv(post_file, nrows=0)
    cols_to_read = [post_header.columns[0]]
    cols_to_read.extend([col for col in obs_names if col in post_header.columns])

    post_oe = pd.read_csv(post_file, usecols=cols_to_read, index_col=0)
    print(f"  Loaded {len(post_oe)} realizations, {len(cols_to_read)-1} observations")

    # Apply filter
    if filtered_ids:
        post_oe = post_oe.loc[[r for r in filtered_ids if r in post_oe.index]]
        print(f"  Filtered to {len(post_oe)} realizations")

    # Initialize arrays for statistics
    stats_names = ['mean', 'median', 'p05', 'p10', 'p25', 'p75', 'p90', 'p95']
    stats_arrays = {name: np.full((nrow, ncol), np.nan) for name in stats_names}

    # Calculate statistics for each cell
    print("\nCalculating statistics...")
    n_cells = 0
    for obs_name, obs_row in arr_h_filtered.iterrows():
        if obs_name not in post_oe.columns:
            continue

        i = int(obs_row['i'])
        j = int(obs_row['j'])

        values = post_oe[obs_name].values

        stats_arrays['mean'][i, j] = np.mean(values)
        stats_arrays['median'][i, j] = np.median(values)
        stats_arrays['p05'][i, j] = np.percentile(values, 5)
        stats_arrays['p10'][i, j] = np.percentile(values, 10)
        stats_arrays['p25'][i, j] = np.percentile(values, 25)
        stats_arrays['p75'][i, j] = np.percentile(values, 75)
        stats_arrays['p90'][i, j] = np.percentile(values, 90)
        stats_arrays['p95'][i, j] = np.percentile(values, 95)

        n_cells += 1

    print(f"  Calculated statistics for {n_cells} cells")

    # Save arrays
    print("\nSaving arrays...")
    for name, arr in stats_arrays.items():
        filename = f'{run_name}_heads_kper{kper}_kstp{kstp}_{name}{file_suffix}.npy'
        filepath = os.path.join(arrays_dir, filename)
        np.save(filepath, arr)
        print(f"  Saved: {filepath}")

    # Also save as CSV for easier viewing
    for name, arr in stats_arrays.items():
        filename = f'{run_name}_heads_kper{kper}_kstp{kstp}_{name}{file_suffix}.csv'
        filepath = os.path.join(arrays_dir, filename)
        pd.DataFrame(arr).to_csv(filepath, index=False, header=False)

    # Create visualization plots
    print("\nCreating plots...")

    # Plot 1: Array maps for each statistic
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes_flat = axes.flatten()

    for idx, (name, arr) in enumerate(stats_arrays.items()):
        ax = axes_flat[idx]

        # Mask NaN values
        masked_arr = np.ma.masked_invalid(arr)

        im = ax.imshow(masked_arr, cmap='viridis', aspect='equal', vmin=2, vmax=20)
        ax.set_title(name, fontsize=10, fontweight='bold')
        ax.set_xlabel('Column', fontsize=8)
        ax.set_ylabel('Row', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.tick_params(labelsize=7)

    plt.suptitle(f'Head Statistics - kper={kper}, kstp={kstp}', fontsize=12, fontweight='bold')
    plt.tight_layout()

    plot_path = os.path.join(output_dir, f'{run_name}_heads_kper{kper}_kstp{kstp}_arrays{file_suffix}.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_path}")

    # Plot 2: Histograms of each statistic
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes_flat = axes.flatten()

    for idx, (name, arr) in enumerate(stats_arrays.items()):
        ax = axes_flat[idx]

        # Get non-NaN values
        valid_vals = arr[~np.isnan(arr)].flatten()

        if len(valid_vals) > 0:
            ax.hist(valid_vals, bins=50, color='#1f77b4', edgecolor='black', linewidth=0.3)

            # Add vertical lines for mean and median
            arr_mean = np.mean(valid_vals)
            arr_median = np.median(valid_vals)
            ax.axvline(arr_mean, color='red', linestyle='--', linewidth=1.5, label=f'Mean: {arr_mean:.2f}')
            ax.axvline(arr_median, color='orange', linestyle='-', linewidth=1.5, label=f'Median: {arr_median:.2f}')
            ax.legend(fontsize=6)

        ax.set_title(name, fontsize=10, fontweight='bold')
        ax.set_xlabel('Head (m)', fontsize=8)
        ax.set_ylabel('Frequency', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(True, alpha=0.3)

    plt.suptitle(f'Head Statistics Histograms - kper={kper}, kstp={kstp}', fontsize=12, fontweight='bold')
    plt.tight_layout()

    hist_path = os.path.join(output_dir, f'{run_name}_heads_kper{kper}_kstp{kstp}_histograms{file_suffix}.png')
    plt.savefig(hist_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {hist_path}")

    # Plot 3: Uncertainty range (p95 - p05)
    fig, ax = plt.subplots(1, 1, figsize=(10, 8))

    uncertainty = stats_arrays['p95'] - stats_arrays['p05']
    masked_unc = np.ma.masked_invalid(uncertainty)

    im = ax.imshow(masked_unc, cmap='YlOrRd', aspect='equal')
    ax.set_title(f'Head Uncertainty (P95 - P05) - kper={kper}, kstp={kstp}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Column', fontsize=10)
    ax.set_ylabel('Row', fontsize=10)
    plt.colorbar(im, ax=ax, label='Uncertainty (m)')

    unc_path = os.path.join(output_dir, f'{run_name}_heads_kper{kper}_kstp{kstp}_uncertainty{file_suffix}.png')
    plt.savefig(unc_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {unc_path}")

    # Save summary statistics
    summary_data = {
        'statistic': stats_names,
        'min': [np.nanmin(stats_arrays[name]) for name in stats_names],
        'max': [np.nanmax(stats_arrays[name]) for name in stats_names],
        'mean': [np.nanmean(stats_arrays[name]) for name in stats_names],
        'std': [np.nanstd(stats_arrays[name]) for name in stats_names],
    }
    summary_df = pd.DataFrame(summary_data)
    summary_path = os.path.join(output_dir, f'{run_name}_heads_kper{kper}_kstp{kstp}_summary{file_suffix}.csv')
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

    create_heads_stats(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        kper=args.kper,
        kstp=args.kstp,
        filter_file=args.filter_file,
        suffix=args.suffix
    )
