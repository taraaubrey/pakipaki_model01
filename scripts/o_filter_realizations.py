"""
o_filter_realizations.py - Filter ensemble realizations based on various criteria

Creates a filtered list of realization indices for subsequent analysis.

Filter options (any combination):
- Iteration selection (default: last)
- Phi percentile threshold
- Head constraints (max mean head threshold)
- Confined flux constraints (no fluxes > 0)

Usage:
    python scripts/o_filter_realizations.py run_name [options]

Examples:
    python scripts/o_filter_realizations.py local_run34 --iteration 19
    python scripts/o_filter_realizations.py local_run34 -i 19 --max-head 20 --phi-percentile 90
    python scripts/o_filter_realizations.py local_run34 -i 19 --max-head 20 --no-positive-confined-flux
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import pyemu

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import default model name from setup
try:
    from setup import MODEL_NAME as DEFAULT_MODEL_NAME
except ImportError:
    DEFAULT_MODEL_NAME = None


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Filter ensemble realizations based on various criteria',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default=DEFAULT_MODEL_NAME,
        help=f'Model name for file prefixes (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )

    parser.add_argument(
        '--run-name', '-r',
        type=str,
        default=None,
        help='Run name for directory path (default: same as model_name)'
    )

    parser.add_argument(
        '--iteration', '-i',
        type=int,
        default=None,
        help='Iteration to analyze (default: last available)'
    )

    parser.add_argument(
        '--phi-percentile',
        type=float,
        nargs=2,
        metavar=('MIN', 'MAX'),
        default=None,
        help='Keep realizations with phi between these percentiles (e.g., 5 95 keeps middle 90%%)'
    )

    parser.add_argument(
        '--max-head',
        type=float,
        default=None,
        help='Maximum mean head threshold - filter out realizations with mean head above this value (e.g., 20)'
    )

    parser.add_argument(
        '--no-positive-confined-flux',
        action='store_true',
        help='Filter out realizations with any positive confined flux (flux > 0)'
    )

    parser.add_argument(
        '--output', '-o',
        type=str,
        default=None,
        help='Output file path for filtered indices (default: models/{run}/filtered_realizations.csv)'
    )

    parser.add_argument(
        '--suffix',
        type=str,
        default=None,
        help='Suffix for output file (e.g., "_filtered" -> filtered_realizations_filtered.csv)'
    )

    return parser.parse_args()


def get_last_iteration(master_dir, model_name):
    """Auto-detect the last available iteration."""
    obs_files = [f for f in os.listdir(master_dir) if f.endswith('.obs.csv')]
    iters = []
    for f in obs_files:
        try:
            iter_num = int(f.split('.')[1])
            iters.append(iter_num)
        except (ValueError, IndexError):
            pass
    return max(iters) if iters else 0


def filter_by_phi(oe, pst, master_dir, model_name, iteration, percentile_range):
    """
    Filter ensemble by phi value percentile range.

    Parameters:
    -----------
    percentile_range : list
        [min_percentile, max_percentile] - keep realizations between these percentiles

    Returns mask of realizations to keep.
    """
    # Load phi values for this iteration
    phi_file = os.path.join(master_dir, f'{model_name}.{iteration}.phi.actual.csv')

    if not os.path.exists(phi_file):
        # Try alternative phi file
        phi_file = os.path.join(master_dir, f'{model_name}.phi.actual.csv')

    if os.path.exists(phi_file):
        phi_df = pd.read_csv(phi_file, index_col=0)
        # Get phi values for the iteration
        if iteration in phi_df.columns:
            phi_values = phi_df[iteration]
        elif str(iteration) in phi_df.columns:
            phi_values = phi_df[str(iteration)]
        else:
            # Assume it's a single iteration file with 'total' column
            if 'total' in phi_df.columns:
                phi_values = phi_df['total']
            else:
                # Use first column
                phi_values = phi_df.iloc[:, 0]

        # Calculate percentile thresholds
        min_pct, max_pct = percentile_range
        threshold_low = np.percentile(phi_values.dropna(), min_pct)
        threshold_high = np.percentile(phi_values.dropna(), max_pct)
        mask = (phi_values >= threshold_low) & (phi_values <= threshold_high)

        # Align with observation ensemble index
        mask = mask.reindex(oe.index, fill_value=False)

        print(f"  Filter (phi between {min_pct}th-{max_pct}th percentile = {threshold_low:.2f}-{threshold_high:.2f}): {mask.sum()} realizations pass")
        return mask
    else:
        print(f"  Warning: Phi file not found, skipping phi filter")
        return pd.Series([True] * len(oe), index=oe.index)


def filter_by_max_head(oe, obs_data, max_head):
    """
    Filter ensemble by maximum mean head threshold.
    Removes realizations where mean head observations exceed threshold.

    Returns mask of realizations to keep.
    """
    # Get head observation names (arr-h type observations)
    head_obs = obs_data[obs_data['oname'].str.contains('arr-h|head', case=False, na=False)].index.tolist()

    if len(head_obs) == 0:
        # Try alternative: look for observations that might be heads
        head_obs = obs_data[obs_data.index.str.contains('head', case=False)].index.tolist()

    if len(head_obs) > 0:
        head_cols = [col for col in oe.columns if col in head_obs]
        if len(head_cols) > 0:
            # Calculate mean head for each realization
            mean_head = oe[head_cols].mean(axis=1)
            mask = mean_head <= max_head
            print(f"  Filter (mean head <= {max_head}): {mask.sum()} realizations pass")
            print(f"    Mean head range: {mean_head.min():.2f} to {mean_head.max():.2f}")
            return mask
        else:
            print(f"  Warning: No head observations found in ensemble, skipping head filter")
            return pd.Series([True] * len(oe), index=oe.index)
    else:
        print(f"  Warning: No head observations found, skipping head filter")
        return pd.Series([True] * len(oe), index=oe.index)


def filter_by_confined_flux(oe, obs_data):
    """
    Filter ensemble by confined flux constraint.
    Removes realizations with any positive confined flux.

    Returns mask of realizations to keep.
    """
    # Get confined budget/flux observation names
    conf_obs = obs_data[
        obs_data['oname'].str.contains('budget.*confined|confined.*flux|conf.*budget', case=False, na=False) |
        obs_data.index.str.contains('conf.*budget|budget.*conf', case=False)
    ].index.tolist()

    if len(conf_obs) == 0:
        # Try alternative patterns
        conf_obs = obs_data[obs_data.index.str.contains('ghb_conf', case=False)].index.tolist()

    if len(conf_obs) > 0:
        conf_cols = [col for col in oe.columns if col in conf_obs]
        if len(conf_cols) > 0:
            # Filter out realizations with ANY positive confined flux
            mask = (oe[conf_cols] <= 0).all(axis=1)
            print(f"  Filter (confined flux <= 0): {mask.sum()} realizations pass")
            return mask
        else:
            print(f"  Warning: No confined flux observations found in ensemble, skipping flux filter")
            return pd.Series([True] * len(oe), index=oe.index)
    else:
        print(f"  Warning: No confined flux observations found, skipping flux filter")
        return pd.Series([True] * len(oe), index=oe.index)


def filter_realizations(run_name, model_name, iteration=None, phi_percentile=None,
                        max_head=None, no_positive_confined_flux=False, output_file=None, suffix=None):
    """
    Filter realizations based on specified criteria and save filtered indices.

    Parameters:
    -----------
    run_name : str
        Name of the run directory
    model_name : str
        Name of the model
    iteration : int, optional
        Iteration to analyze (default: last)
    phi_percentile : float, optional
        Keep realizations below this phi percentile
    max_head : float, optional
        Maximum mean head threshold
    no_positive_confined_flux : bool
        If True, filter out realizations with positive confined flux
    output_file : str, optional
        Path to save filtered indices
    suffix : str, optional
        Suffix for output file (e.g., "_filtered")

    Returns:
    --------
    list : Filtered realization indices
    """

    print(f"\n{'='*80}")
    print(f"FILTER REALIZATIONS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')

    if not os.path.exists(master_dir):
        print(f"Error: Directory not found: {master_dir}")
        return []

    # Auto-detect iteration if not specified
    if iteration is None:
        iteration = get_last_iteration(master_dir, model_name)
        print(f"Using last available iteration: {iteration}")
    else:
        print(f"Iteration: {iteration}")

    # Load PST and observation data
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Load observation ensemble
    oe_file = os.path.join(master_dir, f'{model_name}.{iteration}.obs.csv')
    if not os.path.exists(oe_file):
        print(f"Error: Observation ensemble not found: {oe_file}")
        return []

    oe = pd.read_csv(oe_file, index_col=0)
    print(f"\nInitial ensemble size: {len(oe)} realizations")

    # Track filters applied
    filters_applied = []

    # Initialize combined mask
    mask_combined = pd.Series([True] * len(oe), index=oe.index)

    # Apply filters
    print("\nApplying filters...")

    # Filter 1: Phi percentile range
    if phi_percentile is not None:
        mask_phi = filter_by_phi(oe, pst, master_dir, model_name, iteration, phi_percentile)
        mask_combined = mask_combined & mask_phi
        filters_applied.append(f'phi_percentile_{phi_percentile[0]}_{phi_percentile[1]}')

    # Filter 2: Max head
    if max_head is not None:
        mask_head = filter_by_max_head(oe, obs_data, max_head)
        mask_combined = mask_combined & mask_head
        filters_applied.append(f'max_head_{max_head}')

    # Filter 3: Confined flux
    if no_positive_confined_flux:
        mask_flux = filter_by_confined_flux(oe, obs_data)
        mask_combined = mask_combined & mask_flux
        filters_applied.append('no_positive_confined_flux')

    # Get filtered indices
    filtered_indices = oe.index[mask_combined].tolist()

    print(f"\n{'='*60}")
    print(f"FILTERING RESULTS")
    print(f"{'='*60}")
    print(f"  Initial: {len(oe)} realizations")
    print(f"  Filtered: {len(filtered_indices)} realizations")
    print(f"  Removed: {len(oe) - len(filtered_indices)} realizations")
    print(f"  Retention: {len(filtered_indices)/len(oe)*100:.1f}%")

    if len(filters_applied) == 0:
        print(f"\n  Note: No filters applied, all realizations kept")
        filters_applied.append('none')
    else:
        print(f"\n  Filters applied: {', '.join(filters_applied)}")

    # Save filtered indices
    file_suffix = suffix if suffix else ''
    if output_file is None:
        output_file = os.path.join('models', run_name, f'filtered_realizations{file_suffix}.csv')

    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Save as DataFrame with metadata
    output_df = pd.DataFrame({
        'realization': filtered_indices
    })

    # Add metadata as comments in file header
    with open(output_file, 'w') as f:
        f.write(f"# Filtered realizations for {run_name}\n")
        f.write(f"# Iteration: {iteration}\n")
        f.write(f"# Filters: {', '.join(filters_applied)}\n")
        f.write(f"# Initial count: {len(oe)}\n")
        f.write(f"# Filtered count: {len(filtered_indices)}\n")
        output_df.to_csv(f, index=False)

    print(f"\n  Saved filtered indices to: {output_file}")

    return filtered_indices


def load_filtered_realizations(run_name):
    """
    Load previously filtered realization indices.

    Parameters:
    -----------
    run_name : str
        Name of the run directory

    Returns:
    --------
    list : Filtered realization indices, or None if file not found
    """
    filter_file = os.path.join('models', run_name, 'filtered_realizations.csv')

    if not os.path.exists(filter_file):
        return None

    # Read file, skipping comment lines
    df = pd.read_csv(filter_file, comment='#')
    return df['realization'].tolist()


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    filtered = filter_realizations(
        run_name=run_name,
        model_name=model_name,
        iteration=args.iteration,
        phi_percentile=args.phi_percentile,
        max_head=args.max_head,
        no_positive_confined_flux=args.no_positive_confined_flux,
        output_file=args.output,
        suffix=args.suffix
    )

    if len(filtered) > 0:
        print(f"\nFiltered realization indices: {filtered[:10]}{'...' if len(filtered) > 10 else ''}")
