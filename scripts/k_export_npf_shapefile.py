"""
Export NPF (hydraulic conductivity) statistics to shapefile.

This script:
1. Loads a parameter ensemble (e.g., subset_pe_phi10-90.csv)
2. Extracts NPF parameters for kper1 kstp1
3. Calculates statistics: mean, CV, min, max, std_log
4. Joins to grid shapefile
5. Saves as new shapefile in spatial folder

Usage:
    python scripts/k_export_npf_shapefile.py local_run32 figures/subset_pe_phi10-90.csv
    python scripts/k_export_npf_shapefile.py local_run32 --pe-file figures/subset_pe_phi10-90.csv
    python scripts/k_export_npf_shapefile.py local_run32 --last-iter 16  # Use full posterior
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
import geopandas as gpd
from setup import MODEL_NAME


def calculate_stats(values):
    """Calculate statistics for a parameter across realizations."""
    values = np.array(values)

    # Filter out invalid values
    valid_values = values[~np.isnan(values)]

    if len(valid_values) == 0:
        return {
            'mean': np.nan,
            'cv': np.nan,
            'min': np.nan,
            'max': np.nan,
            'std_log': np.nan
        }

    mean_val = np.mean(valid_values)
    std_val = np.std(valid_values)

    # Calculate CV (coefficient of variation)
    cv = (std_val / mean_val * 100) if mean_val != 0 else np.nan

    # Calculate std of log-transformed values
    try:
        log_values = np.log10(valid_values[valid_values > 0])
        std_log = np.std(log_values) if len(log_values) > 0 else np.nan
    except:
        std_log = np.nan

    return {
        'mean': mean_val,
        'cv': cv,
        'min': np.min(valid_values),
        'max': np.max(valid_values),
        'std_log': std_log
    }


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Export NPF statistics to shapefile',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default=MODEL_NAME,
        help=f'Model/run name (default: {MODEL_NAME} from setup.py)'
    )

    parser.add_argument(
        '--pe-file',
        type=str,
        default=None,
        help='Path to parameter ensemble CSV file (e.g., figures/subset_pe_phi10-90.csv)'
    )

    parser.add_argument(
        '--last-iter',
        type=int,
        default=None,
        help='Use full posterior ensemble from this iteration'
    )

    parser.add_argument(
        '--kper',
        type=int,
        default=1,
        help='Stress period (default: 1)'
    )

    parser.add_argument(
        '--kstp',
        type=int,
        default=1,
        help='Time step (default: 1)'
    )

    parser.add_argument(
        '--output-name',
        type=str,
        default=None,
        help='Output shapefile name (default: auto-generated)'
    )

    return parser.parse_args()


def main():
    args = parse_args()
    model_name = args.model_name

    print(f"\n{'='*80}")
    print(f"EXPORT NPF STATISTICS TO SHAPEFILE: {model_name}")
    print(f"{'='*80}\n")

    # Determine model directory
    model_ws = os.path.join('models', model_name)
    if not os.path.exists(model_ws):
        print(f"ERROR: Model directory not found: {model_ws}")
        return 1

    spatial_dir = os.path.join(model_ws, 'spatial')
    if not os.path.exists(spatial_dir):
        print(f"ERROR: Spatial directory not found: {spatial_dir}")
        return 1

    # Load grid shapefile
    grid_shp = os.path.join(spatial_dir, f'{model_name}_grid.shp')
    if not os.path.exists(grid_shp):
        print(f"ERROR: Grid shapefile not found: {grid_shp}")
        return 1

    print(f"Loading grid shapefile: {grid_shp}")
    grid = gpd.read_file(grid_shp)
    print(f"Grid loaded: {len(grid)} cells")
    print(f"Grid columns: {grid.columns.tolist()}")

    # Determine parameter ensemble file
    if args.pe_file:
        pe_file = os.path.join(model_ws, args.pe_file) if not os.path.isabs(args.pe_file) else args.pe_file
    elif args.last_iter is not None:
        pe_file = os.path.join(model_ws, 'pest', 'master_ies', f'{model_name}.{args.last_iter}.par.csv')
    else:
        print("ERROR: Must specify either --pe-file or --last-iter")
        return 1

    if not os.path.exists(pe_file):
        print(f"ERROR: Parameter ensemble file not found: {pe_file}")
        return 1

    print(f"\nLoading parameter ensemble: {pe_file}")
    pe = pd.read_csv(pe_file, index_col=0)
    print(f"Parameter ensemble loaded: {pe.shape[0]} realizations, {pe.shape[1]} parameters")

    # Load parameter data to get actual parameter names and spatial indices
    par_data_file = os.path.join(model_ws, 'pest', 'master_ies', f'{model_name}.par_data.csv')
    if not os.path.exists(par_data_file):
        print(f"ERROR: Parameter data file not found: {par_data_file}")
        return 1

    print(f"Loading parameter data: {par_data_file}")
    par_data = pd.read_csv(par_data_file, index_col=0)
    print(f"Parameter data loaded: {len(par_data)} parameters")

    # Filter for NPF parameters only
    # NPF parameters have pname like 'npfklayer1-gr', 'npfklayer1-pp', 'npfklayer1-cn'
    npf_params = par_data[par_data['pname'].str.contains('npfklayer', case=False, na=False)]

    if len(npf_params) == 0:
        print("ERROR: No NPF parameters found in parameter data")
        print("Available pnames:")
        print(f"  {par_data['pname'].unique()}")
        return 1

    print(f"\nFound {len(npf_params)} NPF parameters")
    print(f"Parameter types: {npf_params['pname'].unique()}")

    # Extract parameter values and calculate statistics
    param_data = []

    for parnme, row in npf_params.iterrows():
        # Check if this parameter exists in the ensemble
        if parnme not in pe.columns:
            continue

        i_val = row['i']
        j_val = row['j']

        # Skip if i,j are NaN
        if pd.isna(i_val) or pd.isna(j_val):
            continue

        # Get parameter values across all realizations
        values = pe[parnme].values

        # Calculate statistics
        stats = calculate_stats(values)

        param_data.append({
            'row': int(i_val),
            'col': int(j_val),
            'parnme': parnme,
            'pname': row['pname'],
            **stats
        })

    if len(param_data) == 0:
        print("ERROR: Could not extract any NPF parameter statistics")
        return 1

    print(f"\nSuccessfully extracted statistics for {len(param_data)} NPF parameters")

    # Create DataFrame with statistics
    stats_df = pd.DataFrame(param_data)
    print(f"\nStatistics calculated:")
    print(f"  Mean K range: {stats_df['mean'].min():.2e} to {stats_df['mean'].max():.2e}")
    print(f"  CV range: {stats_df['cv'].min():.2f}% to {stats_df['cv'].max():.2f}%")

    # Join to grid based on row, col indices
    if 'row' not in grid.columns or 'col' not in grid.columns:
        print(f"ERROR: Grid does not have row/col columns")
        print(f"Available columns: {grid.columns.tolist()}")
        return 1

    print(f"\nJoining statistics to grid on (row, col)")

    # Merge
    grid_with_stats = grid.merge(
        stats_df[['row', 'col', 'mean', 'cv', 'min', 'max', 'std_log', 'pname']],
        on=['row', 'col'],
        how='left'
    )

    # Rename columns for clarity
    grid_with_stats = grid_with_stats.rename(columns={
        'mean': 'K_mean',
        'cv': 'K_cv',
        'min': 'K_min',
        'max': 'K_max',
        'std_log': 'K_std_log',
        'pname': 'K_pname'
    })

    # Generate output filename
    if args.output_name:
        output_shp = os.path.join(spatial_dir, args.output_name)
        if not output_shp.endswith('.shp'):
            output_shp += '.shp'
    else:
        # Auto-generate name based on input file
        if args.pe_file:
            pe_basename = os.path.basename(args.pe_file).replace('.csv', '')
        else:
            pe_basename = f'iter{args.last_iter}'
        output_shp = os.path.join(spatial_dir, f'{model_name}_npf_kper{args.kper}_kstp{args.kstp}_{pe_basename}.shp')

    # Save shapefile
    print(f"\nSaving shapefile: {output_shp}")
    grid_with_stats.to_file(output_shp)

    print(f"\n{'='*80}")
    print("EXPORT COMPLETE!")
    print(f"{'='*80}")
    print(f"Output: {output_shp}")
    print(f"Cells with data: {grid_with_stats['K_mean'].notna().sum()} / {len(grid_with_stats)}")
    print(f"\nStatistics summary:")
    print(grid_with_stats[['K_mean', 'K_cv', 'K_min', 'K_max', 'K_std_log']].describe())
    print(f"{'='*80}\n")

    return 0


if __name__ == '__main__':
    sys.exit(main())
