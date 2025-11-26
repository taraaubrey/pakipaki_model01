"""
pp_springflux_vs_phi_dask.py - Plot spring flux predictions vs objective function (Dask version)

Scatter plot comparing spring flux predictions (kper 3, past steady state)
against objective function values for each realization across iterations.

Uses Dask for parallel processing of large observation files.

Usage:
    python scripts/pp_springflux_vs_phi_dask.py run_name [options]

Examples:
    python scripts/pp_springflux_vs_phi_dask.py run37 --run-name run37_62894 --n-workers 120
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemu
from concurrent.futures import ProcessPoolExecutor, as_completed

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
        description='Plot spring flux predictions vs objective function (Dask version)'
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
    parser.add_argument('--iterations', type=str, default='0,19,39',
                       help='Comma-separated iterations to plot (default: 0,19,39)')
    parser.add_argument('--filter-file', type=str, default=None,
                       help='Path to file with filtered realization indices')
    parser.add_argument('--suffix', type=str, default='',
                       help='Suffix for output files')
    parser.add_argument('--n-workers', type=int, default=None,
                       help='Number of Dask workers (default: number of CPUs)')
    parser.add_argument('--use-dask', action='store_true',
                       help='Use Dask distributed for parallel processing')
    return parser.parse_args()


def load_filtered_realizations(filter_file):
    """Load filtered realization indices from file."""
    if filter_file is None or not os.path.exists(filter_file):
        return None
    df = pd.read_csv(filter_file, comment='#')
    return df['realization'].tolist()


def process_iteration(args):
    """Process a single iteration - designed for parallel execution."""
    iter_num, master_dir, model_name, spring_obs, filtered_ids = args

    obs_file = os.path.join(master_dir, f'{model_name}.{iter_num}.obs.csv')
    if not os.path.exists(obs_file):
        return iter_num, None

    # Read only the column we need for efficiency
    try:
        # First get the header to find column index
        header = pd.read_csv(obs_file, nrows=0)
        if spring_obs not in header.columns:
            return iter_num, None

        # Read only index and target column
        cols_to_read = [header.columns[0], spring_obs]
        obs_df = pd.read_csv(obs_file, usecols=cols_to_read, index_col=0, low_memory=False)

        spring_vals = obs_df[spring_obs]

        # Convert index to strings for consistent matching
        result = {
            'index': [str(i) for i in spring_vals.index],
            'values': spring_vals.values
        }

        return iter_num, result

    except Exception as e:
        print(f"  Error processing iteration {iter_num}: {e}")
        return iter_num, None


def analyze_springflux_vs_phi_dask(run_name, model_name, iterations=[0, 19, 39],
                                    filter_file=None, suffix='', n_workers=None, use_dask=False):
    """
    Create scatter plot of spring flux vs objective function using Dask.
    """
    print(f"\n{'='*80}")
    print(f"SPRING FLUX VS PHI ANALYSIS (DASK): {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_springflux_vs_phi')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Determine number of workers
    if n_workers is None:
        n_workers = os.cpu_count()
    print(f"\nUsing {n_workers} workers")

    # Load PST to get observation info
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    print(f"\nLoading PST: {pst_file}")
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Find spring flux observation for kper 3 (past steady state)
    spring_kper3_obs = None

    # Look for GHB-Spring budget observation at kper 3
    for idx, row in obs_data.iterrows():
        oname = str(row.get('oname', ''))
        obgnme = str(row.get('obgnme', ''))
        if 'ghb-sp' in oname.lower() or 'ghb_sp' in obgnme.lower():
            kper = row.get('kper', None)
            if kper is not None and int(kper) == 3:
                spring_kper3_obs = idx
                break

    if spring_kper3_obs is None:
        # Fall back to arr-spq
        arr_spq = obs_data[obs_data['oname'] == 'arr-spq']
        if len(arr_spq) > 0:
            arr_spq_kper3 = arr_spq[arr_spq['kper'].astype(int) == 3]
            if len(arr_spq_kper3) > 0:
                spring_kper3_obs = arr_spq_kper3.index[0]
                print(f"Using arr-spq observation: {spring_kper3_obs}")

    if spring_kper3_obs is None:
        print("Could not find spring flux observation for kper 3")
        return

    print(f"\nSpring flux observation (kper 3): {spring_kper3_obs}")

    # Load phi values
    phi_file = os.path.join(master_dir, f'{model_name}.phi.actual.csv')
    print(f"Loading phi values: {phi_file}")

    if not os.path.exists(phi_file):
        print(f"  Phi file not found: {phi_file}")
        return

    phi_df = pd.read_csv(phi_file, index_col=0)
    print(f"  Phi data shape: {phi_df.shape}")

    # Load filtered realizations
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations")

    # Process iterations in parallel
    print(f"\nLoading observation data for {len(iterations)} iterations...")

    if use_dask:
        # Use Dask distributed
        try:
            from dask.distributed import Client, LocalCluster

            print(f"Starting Dask cluster with {n_workers} workers...")
            cluster = LocalCluster(n_workers=n_workers, threads_per_worker=1)
            client = Client(cluster)
            print(f"  Dashboard: {client.dashboard_link}")

            # Submit tasks
            futures = []
            for iter_num in iterations:
                args = (iter_num, master_dir, model_name, spring_kper3_obs, filtered_ids)
                future = client.submit(process_iteration, args)
                futures.append(future)

            # Gather results
            obs_results = {}
            for future in as_completed(futures):
                iter_num, result = future.result()
                if result is not None:
                    obs_results[iter_num] = result
                    print(f"  Loaded iteration {iter_num}: {len(result['values'])} realizations")

            client.close()
            cluster.close()

        except ImportError:
            print("Dask distributed not available, falling back to ProcessPoolExecutor")
            use_dask = False

    if not use_dask:
        # Use ProcessPoolExecutor for parallel processing
        obs_results = {}

        # Limit workers for I/O bound tasks
        actual_workers = min(n_workers, len(iterations), 8)  # I/O bound, don't need many

        with ProcessPoolExecutor(max_workers=actual_workers) as executor:
            args_list = [(iter_num, master_dir, model_name, spring_kper3_obs, filtered_ids)
                        for iter_num in iterations]

            future_to_iter = {executor.submit(process_iteration, args): args[0]
                             for args in args_list}

            for future in as_completed(future_to_iter):
                iter_num, result = future.result()
                if result is not None:
                    obs_results[iter_num] = result
                    print(f"  Loaded iteration {iter_num}: {len(result['values'])} realizations")

    # Build plot data
    plot_data = {}

    for iter_num in iterations:
        if iter_num not in obs_results:
            continue

        iter_col = str(iter_num)
        if iter_col not in phi_df.columns:
            print(f"  Iteration {iter_num} not in phi data")
            continue

        spring_index = obs_results[iter_num]['index']
        spring_values = obs_results[iter_num]['values']

        phi_vals = phi_df[iter_col]
        phi_index = [str(i) for i in phi_vals.index]

        # Match realizations
        common_reals = set(spring_index) & set(phi_index)

        if filtered_ids:
            filtered_str = [str(r) for r in filtered_ids]
            common_reals = common_reals & set(filtered_str)

        # Build matched arrays
        spring_matched = []
        phi_matched = []
        real_ids = []

        for r in common_reals:
            if r in spring_index and r in phi_index:
                s_idx = spring_index.index(r)
                p_idx = phi_index.index(r)
                spring_matched.append(spring_values[s_idx])
                phi_matched.append(phi_vals.iloc[p_idx])
                real_ids.append(r)

        if len(spring_matched) > 0:
            plot_data[iter_num] = {
                'spring': np.array(spring_matched),
                'phi': np.array(phi_matched),
                'reals': real_ids
            }

    if len(plot_data) == 0:
        print("\nNo data to plot")
        return

    # Create figure
    print("\nCreating plot...")
    n_iters = len(plot_data)
    fig, axes = plt.subplots(1, n_iters, figsize=(7.28, 7.28/3))

    if n_iters == 1:
        axes = [axes]

    colors = {
        0: '#8d8d8d',    # grey - prior
        19: '#2ca02c',   # green - mid iteration
        39: '#ff7f0e',   # orange - final
    }

    panel_labels = ['A', 'B', 'C', 'D', 'E', 'F']

    for idx, (iter_num, data) in enumerate(sorted(plot_data.items())):
        ax = axes[idx]

        spring = data['spring']
        phi = data['phi']

        color = colors.get(iter_num, '#1f77b4')

        ax.scatter(phi, spring, c=color, s=10, alpha=0.6, edgecolors='none')

        # Add correlation
        corr = np.corrcoef(phi, spring)[0, 1]

        # Labels
        if iter_num == 0:
            title = f"{panel_labels[idx]}: Prior (iter 0)"
        else:
            title = f"{panel_labels[idx]}: Iter {iter_num}"

        ax.set_title(title, fontsize=6, fontweight='bold', loc='left')
        ax.set_xlabel('Objective Function (Phi)', fontsize=6)
        ax.set_ylabel('Spring Flux kper 3 (m³/d)', fontsize=6)
        ax.tick_params(labelsize=5)
        ax.grid(True, alpha=0.3)

        # Add statistics
        stats_text = f"n={len(spring)}\nr={corr:.3f}"
        ax.text(0.98, 0.98, stats_text, transform=ax.transAxes, fontsize=5,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()

    # Save plot
    plot_file = os.path.join(output_dir, f'{run_name}_springflux_vs_phi{suffix}.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved plot to: {plot_file}")

    # Save data
    data_rows = []
    for iter_num, data in plot_data.items():
        for i, r in enumerate(data['reals']):
            data_rows.append({
                'iteration': iter_num,
                'realization': r,
                'spring_flux_kper3': data['spring'][i],
                'phi': data['phi'][i]
            })

    data_df = pd.DataFrame(data_rows)
    data_file = os.path.join(output_dir, f'{run_name}_springflux_vs_phi_data{suffix}.csv')
    data_df.to_csv(data_file, index=False)
    print(f"Saved data to: {data_file}")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    # Parse iterations
    iterations = [int(i.strip()) for i in args.iterations.split(',')]

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")
    print(f"Iterations: {iterations}")

    analyze_springflux_vs_phi_dask(
        run_name=run_name,
        model_name=model_name,
        iterations=iterations,
        filter_file=args.filter_file,
        suffix=args.suffix,
        n_workers=args.n_workers,
        use_dask=args.use_dask
    )
