"""
p_cv_vs_spring_flux.py - Analyze parameter CV evolution vs spring flux predictions

Compares change in CV across iterations for parameters against spring flux predictions
for past (kper 1,2) and future (kper 3,4) periods.

Usage:
    python scripts/p_cv_vs_spring_flux.py run_name [--last-iter N]

Example:
    python scripts/p_cv_vs_spring_flux.py local_run34 --last-iter 19
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemu
import glob

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from setup import REINFLATE_ITERS
from h_analyze_params import categorize_parameter


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Analyze parameter CV evolution vs spring flux predictions'
    )
    parser.add_argument('run_name', type=str, help='Name of the run directory')
    parser.add_argument('--model-name', type=str, default=None,
                       help='Model name (default: same as run_name)')
    parser.add_argument('--last-iter', type=int, default=None,
                       help='Last iteration index (default: auto-detect)')
    parser.add_argument('--spring-i', type=int, default=36,
                       help='Spring cell i index (default: 36)')
    parser.add_argument('--spring-j', type=int, default=33,
                       help='Spring cell j index (default: 33)')
    parser.add_argument('--suffix', type=str, default=None,
                       help='Suffix for output files (e.g., "_filtered")')
    return parser.parse_args()


def get_spring_flux_obs(obs_data, spring_i, spring_j):
    """
    Find spring flux observation names at specified location.

    Returns dict with 'past' and 'future' observation lists.
    """
    # Filter arr-spq observations
    arr_spq = obs_data[obs_data['oname'] == 'arr-spq'].copy()
    arr_spq['i'] = pd.to_numeric(arr_spq['i'], errors='coerce')
    arr_spq['j'] = pd.to_numeric(arr_spq['j'], errors='coerce')
    arr_spq['kper'] = pd.to_numeric(arr_spq['kper'], errors='coerce')

    spring_obs = arr_spq[(arr_spq['i'] == spring_i) & (arr_spq['j'] == spring_j)]

    # Separate by stress period (past: 1,2; future: 3,4)
    past_obs = spring_obs[spring_obs['kper'].isin([1, 2])].index.tolist()
    future_obs = spring_obs[spring_obs['kper'].isin([3, 4])].index.tolist()

    return {'past': past_obs, 'future': future_obs}


def load_spring_flux_cv(master_dir, model_name, spring_obs_names, iterations):
    """
    Load spring flux CV for each iteration by reading only the needed columns.

    Returns DataFrame with iterations as index and CV values.
    """
    results = []

    for iter_num in iterations:
        obs_file = os.path.join(master_dir, f'{model_name}.{iter_num}.obs.csv')
        if not os.path.exists(obs_file):
            continue

        # Read only the spring flux columns (more efficient)
        try:
            # First read header to get column indices
            header = pd.read_csv(obs_file, nrows=0)
            cols_to_read = ['real_name'] if 'real_name' in header.columns else [header.columns[0]]
            cols_to_read.extend([col for col in spring_obs_names if col in header.columns])

            if len(cols_to_read) <= 1:
                continue

            # Read only needed columns
            df = pd.read_csv(obs_file, usecols=cols_to_read, index_col=0)

            # Calculate CV for spring flux
            values = df.values.flatten()
            mean_val = np.mean(values)
            std_val = np.std(values)
            cv = (std_val / abs(mean_val) * 100) if mean_val != 0 else np.nan

            results.append({
                'iteration': iter_num,
                'mean': mean_val,
                'std': std_val,
                'cv': cv
            })
        except Exception as e:
            print(f"  Warning: Could not load iteration {iter_num}: {e}")
            continue

    return pd.DataFrame(results)


def load_parameter_cv_evolution(master_dir, model_name, pst, iterations):
    """
    Load parameter CV evolution across iterations.

    Returns DataFrame with CV by category for each iteration.
    """
    param_data = pst.parameter_data.copy()

    # Create pargp to pname mapping
    pargp_to_pname = {}
    for pargp in param_data['pargp'].unique():
        pars = param_data[param_data['pargp'] == pargp]
        if len(pars) > 0 and 'pname' in pars.columns:
            pargp_to_pname[pargp] = pars['pname'].iloc[0]
        else:
            pargp_to_pname[pargp] = pargp

    results = []

    for iter_num in iterations:
        par_file = os.path.join(master_dir, f'{model_name}.{iter_num}.par.csv')
        if not os.path.exists(par_file):
            continue

        par_df = pd.read_csv(par_file, index_col=0, low_memory=False)

        # Calculate CV by category
        category_stats = {}

        for pargp in param_data['pargp'].unique():
            group_params = param_data[param_data['pargp'] == pargp].index.tolist()
            group_params = [p for p in group_params if p in par_df.columns]

            if len(group_params) == 0:
                continue

            values = par_df[group_params].values.flatten()
            mean_val = np.mean(values)
            std_val = np.std(values)
            cv = (std_val / mean_val * 100) if mean_val != 0 else np.nan

            # Get category
            pname = pargp_to_pname.get(pargp, pargp)
            category, subcat = categorize_parameter(pname)

            if category not in category_stats:
                category_stats[category] = []
            category_stats[category].append(cv)

        # Average CV by category
        row = {'iteration': iter_num}
        for cat, cvs in category_stats.items():
            row[f'{cat}_cv'] = np.mean(cvs)

        # Also calculate overall mean CV
        all_cvs = [cv for cvs in category_stats.values() for cv in cvs]
        row['overall_cv'] = np.mean(all_cvs) if all_cvs else np.nan

        results.append(row)

    return pd.DataFrame(results)


def create_cv_vs_spring_plots(run_name, model_name, last_iter=None, spring_i=36, spring_j=33, suffix=None):
    """
    Create plots comparing parameter CV evolution to spring flux predictions.
    """
    # Set suffix for output files
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"CV vs SPRING FLUX ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_param_analysis')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Get spring flux observations
    spring_obs = get_spring_flux_obs(obs_data, spring_i, spring_j)
    print(f"\nSpring flux observations at i={spring_i}, j={spring_j}:")
    print(f"  Past (kper 1,2): {len(spring_obs['past'])} observations")
    print(f"  Future (kper 3,4): {len(spring_obs['future'])} observations")

    if len(spring_obs['past']) == 0 and len(spring_obs['future']) == 0:
        print("Error: No spring flux observations found")
        return

    # Find available iterations
    par_files = glob.glob(os.path.join(master_dir, f'{model_name}.*.par.csv'))
    iterations = sorted(set([
        int(os.path.basename(f).split('.')[1])
        for f in par_files
        if os.path.basename(f).split('.')[1].isdigit()
    ]))

    if last_iter is not None:
        iterations = [i for i in iterations if i <= last_iter]

    print(f"\nAnalyzing iterations: {iterations}")

    # Load parameter CV evolution
    print("\nLoading parameter CV evolution...")
    param_cv = load_parameter_cv_evolution(master_dir, model_name, pst, iterations)

    # Load spring flux CV for past and future
    print("Loading spring flux CV (past)...")
    spring_past_cv = load_spring_flux_cv(master_dir, model_name, spring_obs['past'], iterations)

    print("Loading spring flux CV (future)...")
    spring_future_cv = load_spring_flux_cv(master_dir, model_name, spring_obs['future'], iterations)

    # Save extracted data for quick reloading
    print("\nSaving extracted data...")
    data_dir = os.path.join(output_dir, 'data')
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    # Save parameter CV data
    param_cv_file = os.path.join(data_dir, f'{run_name}_param_cv{file_suffix}.csv')
    param_cv.to_csv(param_cv_file, index=False)
    print(f"  Saved: {param_cv_file}")

    # Save spring flux CV data
    spring_past_cv_file = os.path.join(data_dir, f'{run_name}_spring_past_cv{file_suffix}.csv')
    spring_past_cv.to_csv(spring_past_cv_file, index=False)
    print(f"  Saved: {spring_past_cv_file}")

    spring_future_cv_file = os.path.join(data_dir, f'{run_name}_spring_future_cv{file_suffix}.csv')
    spring_future_cv.to_csv(spring_future_cv_file, index=False)
    print(f"  Saved: {spring_future_cv_file}")

    # Create plots
    print("\nCreating plots...")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Colors for categories
    cat_colors = {
        'Hydraulic Conductivity': '#1f77b4',
        'Storage': '#2ca02c',
        'Recharge': '#9467bd',
        'GHB Conductance': '#ff7f0e',
        'overall': '#333333'
    }

    # Plot A: Past spring flux
    ax1 = axes[0]
    ax1_twin = ax1.twinx()

    # Parameter CV lines
    for col in param_cv.columns:
        if col == 'iteration' or not col.endswith('_cv'):
            continue
        cat_name = col.replace('_cv', '')
        color = cat_colors.get(cat_name, '#7f7f7f')
        if cat_name == 'overall':
            ax1.plot(param_cv['iteration'], param_cv[col], '-',
                    color=color, linewidth=2, label=cat_name, alpha=0.8)
        else:
            ax1.plot(param_cv['iteration'], param_cv[col], '--',
                    color=color, linewidth=1.5, label=cat_name, alpha=0.7)

    # Spring flux CV
    if len(spring_past_cv) > 0:
        ax1_twin.plot(spring_past_cv['iteration'], spring_past_cv['cv'],
                     'o-', color='red', linewidth=2, markersize=4, label='Spring flux CV')

    ax1.set_xlabel('Iteration', fontsize=10)
    ax1.set_ylabel('Parameter CV (%)', fontsize=10)
    ax1_twin.set_ylabel('Spring Flux CV (%)', fontsize=10, color='red')
    ax1_twin.tick_params(axis='y', labelcolor='red')
    ax1.set_title('A. Past Spring Flux (kper 1,2)', fontsize=11, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper left', fontsize=8)
    ax1_twin.legend(loc='upper right', fontsize=8)

    # Plot B: Future spring flux
    ax2 = axes[1]
    ax2_twin = ax2.twinx()

    # Parameter CV lines
    for col in param_cv.columns:
        if col == 'iteration' or not col.endswith('_cv'):
            continue
        cat_name = col.replace('_cv', '')
        color = cat_colors.get(cat_name, '#7f7f7f')
        if cat_name == 'overall':
            ax2.plot(param_cv['iteration'], param_cv[col], '-',
                    color=color, linewidth=2, label=cat_name, alpha=0.8)
        else:
            ax2.plot(param_cv['iteration'], param_cv[col], '--',
                    color=color, linewidth=1.5, label=cat_name, alpha=0.7)

    # Spring flux CV
    if len(spring_future_cv) > 0:
        ax2_twin.plot(spring_future_cv['iteration'], spring_future_cv['cv'],
                     'o-', color='red', linewidth=2, markersize=4, label='Spring flux CV')

    ax2.set_xlabel('Iteration', fontsize=10)
    ax2.set_ylabel('Parameter CV (%)', fontsize=10)
    ax2_twin.set_ylabel('Spring Flux CV (%)', fontsize=10, color='red')
    ax2_twin.tick_params(axis='y', labelcolor='red')
    ax2.set_title('B. Future Spring Flux (kper 3,4)', fontsize=11, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper left', fontsize=8)
    ax2_twin.legend(loc='upper right', fontsize=8)

    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, f'{run_name}_cv_vs_spring_flux{file_suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved plot to: {output_path}")


if __name__ == '__main__':
    args = parse_args()

    run_name = args.run_name
    model_name = args.model_name if args.model_name else run_name

    create_cv_vs_spring_plots(
        run_name=run_name,
        model_name=model_name,
        last_iter=args.last_iter,
        spring_i=args.spring_i,
        spring_j=args.spring_j,
        suffix=args.suffix
    )
