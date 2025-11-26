"""
pp_plot_arrays.py - Plot saved array statistics from pp_heads_stats or pp_array_stats

Loads and plots the saved .npy arrays with custom colorbar limits.

Usage:
    python scripts/pp_plot_arrays.py run_name [options]

Examples:
    python scripts/pp_plot_arrays.py local_run34 --kper 1 --kstp 1 --vmin 0 --vmax 20
    python scripts/pp_plot_arrays.py local_run34 --oname arr-spq --kper 1 --kstp 1 --vmin -50 --vmax 0
    python scripts/pp_plot_arrays.py local_run34 --kper 2 --kstp 52 --suffix _filtered --vmin 2 --vmax 15
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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
        description='Plot saved array statistics from pp_heads_stats or pp_array_stats'
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
    parser.add_argument('--oname', type=str, default='arr-h',
                       help='Observation name type (default: arr-h). Options: arr-h, arr-spq, arr-awq')
    parser.add_argument('--kper', type=int, default=1,
                       help='Stress period (default: 1)')
    parser.add_argument('--kstp', type=int, default=1,
                       help='Time step (default: 1)')
    parser.add_argument('--suffix', type=str, default='',
                       help='Suffix for input/output files (e.g., "_filtered")')
    parser.add_argument('--vmin', type=float, default=None,
                       help='Minimum value for colorbar')
    parser.add_argument('--vmax', type=float, default=None,
                       help='Maximum value for colorbar')
    parser.add_argument('--cmap', type=str, default='viridis',
                       help='Colormap name (default: viridis)')
    parser.add_argument('--stats', type=str, nargs='+',
                       default=['mean', 'median', 'p05', 'p10', 'p25', 'p75', 'p90', 'p95'],
                       help='Statistics to plot (default: all)')
    parser.add_argument('--output-suffix', type=str, default='_replot',
                       help='Suffix for output plot file (default: _replot)')
    return parser.parse_args()


def plot_arrays(run_name, model_name, oname='arr-h', kper=1, kstp=1, suffix='',
                vmin=None, vmax=None, cmap='viridis', stats=None, output_suffix='_replot'):
    """
    Load and plot saved array statistics.
    """
    # Clean oname for folder/file naming
    oname_clean = oname.replace('-', '_')

    print(f"\n{'='*80}")
    print(f"PLOT ARRAYS: {run_name}")
    print(f"oname={oname}, kper={kper}, kstp={kstp}")
    print(f"{'='*80}")

    # Setup paths
    if oname == 'arr-h':
        # Check both pp_heads and pp_arr_h for backwards compatibility
        arrays_dir = os.path.join('models', run_name, 'pp_heads', 'arrays')
        if not os.path.exists(arrays_dir):
            arrays_dir = os.path.join('models', run_name, f'pp_{oname_clean}', 'arrays')
        output_dir = os.path.join('models', run_name, 'pp_heads')
    else:
        arrays_dir = os.path.join('models', run_name, f'pp_{oname_clean}', 'arrays')
        output_dir = os.path.join('models', run_name, f'pp_{oname_clean}')

    if not os.path.exists(arrays_dir):
        print(f"Error: Arrays directory not found: {arrays_dir}")
        return

    # Default statistics to plot
    if stats is None:
        stats = ['mean', 'median', 'p05', 'p10', 'p25', 'p75', 'p90', 'p95']

    # Load arrays
    print(f"\nLoading arrays from: {arrays_dir}")
    stats_arrays = {}

    for stat_name in stats:
        # Try different filename patterns
        patterns = [
            f'{run_name}_{oname_clean}_kper{kper}_kstp{kstp}_{stat_name}{suffix}.npy',
            f'{run_name}_heads_kper{kper}_kstp{kstp}_{stat_name}{suffix}.npy',
        ]

        loaded = False
        for pattern in patterns:
            filepath = os.path.join(arrays_dir, pattern)
            if os.path.exists(filepath):
                stats_arrays[stat_name] = np.load(filepath)
                print(f"  Loaded: {pattern}")
                loaded = True
                break

        if not loaded:
            print(f"  Warning: Could not find array for {stat_name}")

    if len(stats_arrays) == 0:
        print("Error: No arrays found to plot")
        print(f"Looking for files like: {run_name}_{oname_clean}_kper{kper}_kstp{kstp}_mean{suffix}.npy")
        return

    # Determine colorbar limits if not specified
    if vmin is None:
        vmin = np.nanmin([np.nanmin(arr) for arr in stats_arrays.values()])
        print(f"  Auto vmin: {vmin:.2f}")
    if vmax is None:
        vmax = np.nanmax([np.nanmax(arr) for arr in stats_arrays.values()])
        print(f"  Auto vmax: {vmax:.2f}")

    # Create plot
    print(f"\nCreating plot with vmin={vmin}, vmax={vmax}, cmap={cmap}...")

    n_stats = len(stats_arrays)
    n_cols = min(4, n_stats)
    n_rows = int(np.ceil(n_stats / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))

    if n_stats == 1:
        axes = np.array([axes])
    axes_flat = axes.flatten() if hasattr(axes, 'flatten') else [axes]

    for idx, (name, arr) in enumerate(stats_arrays.items()):
        if idx >= len(axes_flat):
            break

        ax = axes_flat[idx]

        # Mask NaN values
        masked_arr = np.ma.masked_invalid(arr)

        im = ax.imshow(masked_arr, cmap=cmap, aspect='equal', vmin=vmin, vmax=vmax)
        ax.set_title(name, fontsize=10, fontweight='bold')
        ax.set_xlabel('Column', fontsize=8)
        ax.set_ylabel('Row', fontsize=8)
        plt.colorbar(im, ax=ax, shrink=0.8)
        ax.tick_params(labelsize=7)

    # Hide unused subplots
    for idx in range(len(stats_arrays), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.suptitle(f'{oname} Statistics - kper={kper}, kstp={kstp}\nvmin={vmin}, vmax={vmax}',
                fontsize=12, fontweight='bold')
    plt.tight_layout()

    # Save plot
    plot_filename = f'{run_name}_{oname_clean}_kper{kper}_kstp{kstp}_arrays{suffix}{output_suffix}.png'
    plot_path = os.path.join(output_dir, plot_filename)
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved: {plot_path}")

    # Also create uncertainty plot if we have p95 and p05
    if 'p95' in stats_arrays and 'p05' in stats_arrays:
        fig, ax = plt.subplots(1, 1, figsize=(10, 8))

        uncertainty = stats_arrays['p95'] - stats_arrays['p05']
        masked_unc = np.ma.masked_invalid(uncertainty)

        im = ax.imshow(masked_unc, cmap='YlOrRd', aspect='equal')
        ax.set_title(f'{oname} Uncertainty (P95 - P05) - kper={kper}, kstp={kstp}',
                    fontsize=12, fontweight='bold')
        ax.set_xlabel('Column', fontsize=10)
        ax.set_ylabel('Row', fontsize=10)
        plt.colorbar(im, ax=ax, label='Uncertainty')

        unc_filename = f'{run_name}_{oname_clean}_kper{kper}_kstp{kstp}_uncertainty{suffix}{output_suffix}.png'
        unc_path = os.path.join(output_dir, unc_filename)
        plt.savefig(unc_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"Saved: {unc_path}")

    print(f"\nDone!")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    plot_arrays(
        run_name=run_name,
        model_name=model_name,
        oname=args.oname,
        kper=args.kper,
        kstp=args.kstp,
        suffix=args.suffix,
        vmin=args.vmin,
        vmax=args.vmax,
        cmap=args.cmap,
        stats=args.stats,
        output_suffix=args.output_suffix
    )
