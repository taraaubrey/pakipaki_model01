"""
pp_1to1_by_obsgroup.py - Plot 1:1 observed vs simulated for each observation group

Creates 1:1 plots for each observation group (oname) and unique usecol combination.
Excludes zero-weight observations and constraint observations (less_than/greater_than).

Truth data sources:
- pk4-heads, pk4-diff, pk4-spring-diff, pk4-aw-diff: models/{MODEL_NAME}/truth/output.sample_heads.truth.csv
- aw-q: truth/output.GHB_AW_fluxes.truth.csv (column AWq)
- spring-q: truth/output.GHB_SP_fluxes.truth.csv (column SPq)
- budget-sw, budget-rch: truth/output.budget.truth.csv (columns sw, rch)
- pk4-heads: truth/output.sample_heads.truth.csv
- recession: truth/output.sample_recession_rates.truth.csv
- top-heads: truth/output.heads.truth.dat (excluded from plots)

Usage:
    python scripts/pp_1to1_by_obsgroup.py run_name [options]

Examples:
    python scripts/pp_1to1_by_obsgroup.py run36 --run-name run36_62859_noGHBheads --post-iter 39
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
        description='Plot 1:1 observed vs simulated for each observation group'
    )
    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default='run37',
        help=f'Model name for file prefixes (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )
    parser.add_argument('--run-name', '-r', type=str, default='run37_62894',
                       help='Run name for directory path (default: same as model_name)')
    parser.add_argument('--post-iter', type=int, default=19,
                       help='Posterior iteration (default: 39)')
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


def load_truth_data(run_name):
    """
    Load all truth data files.

    Args:
        run_name: Run name for locating model-specific truth files

    Returns dict of DataFrames: {filename: DataFrame}
    """
    truth_dir = 'truth'
    model_truth_dir = os.path.join('models', run_name, 'truth')

    truth_files = {
        # Model-specific truth files
        'pk4-sample-heads': os.path.join(model_truth_dir, 'output.sample_heads.truth.csv'),

        # Global truth files
        'pk4-heads': os.path.join(truth_dir, 'output.sample_heads.truth.csv'),
        'budget': os.path.join(truth_dir, 'output.budget.truth.csv'),
        'awq': os.path.join(truth_dir, 'output.GHB_AW_fluxes.truth.csv'),
        'spq': os.path.join(truth_dir, 'output.GHB_SP_fluxes.truth.csv'),
        'recession': os.path.join(truth_dir, 'output.sample_recession_rates.truth.csv'),
        'top-heads': os.path.join(truth_dir, 'output.heads.truth.dat'),
    }

    truth_data = {}
    for key, filepath in truth_files.items():
        if os.path.exists(filepath):
            # Handle different file formats
            if filepath.endswith('.dat'):
                # Load .dat files as space-delimited
                truth_data[key] = pd.read_csv(filepath, sep=r'\s+')
            else:
                truth_data[key] = pd.read_csv(filepath)
            print(f"  Loaded {key}: {filepath} ({len(truth_data[key])} rows)")
        else:
            print(f"  Warning: {key} truth file not found: {filepath}")
            truth_data[key] = None

    return truth_data


def get_truth_value(obs_row, truth_data):
    """
    Get the truth value for an observation based on its metadata.

    Args:
        obs_row: Row from obs_data with observation metadata
        truth_data: Dict of truth DataFrames

    Returns:
        Truth value or None if not found
    """
    oname = str(obs_row.get('oname', '')).lower()
    usecol = str(obs_row.get('usecol', ''))
    kper = obs_row.get('kper')
    kstp = obs_row.get('kstp', 1)

    try:
        # Convert kper/kstp to int
        kper = int(kper) if pd.notna(kper) else 1
        kstp = int(kstp) if pd.notna(kstp) else 1

        # pk4-heads - use pk4 column from model-specific truth file
        if oname == 'pk4-heads' and truth_data.get('pk4-sample-heads') is not None:
            df = truth_data['pk4-sample-heads']
            # Map kper/kstp to time index (same logic as pk4-heads)
            if kper == 1:
                time_idx = 1
            elif kper == 2:
                time_idx = 1 + kstp
            elif kper == 3:
                time_idx = 54
            elif kper == 4:
                time_idx = 54 + kstp
            else:
                return None

            # Find row with this time and use pk4 column
            if 'time' in df.columns and 'pk4' in df.columns:
                row = df[df['time'] == time_idx]
                if len(row) > 0:
                    return row['pk4'].values[0]
            return None

        # ts-pk4 observations - check usecol for column name
        elif oname == 'ts-pk4' and truth_data.get('pk4-sample-heads') is not None:
            df = truth_data['pk4-sample-heads']
            if kper == 1:
                time_idx = 1
            elif kper == 2:
                time_idx = 1 + kstp
            elif kper == 3:
                time_idx = 54
            elif kper == 4:
                time_idx = 54 + kstp
            else:
                return None

            # Use usecol as the column name
            if 'time' in df.columns and usecol in df.columns:
                row = df[df['time'] == time_idx]
                if len(row) > 0:
                    return row[usecol].values[0]
            return None

        # pk4-diff - use pk4-diff column (backward compatibility)
        elif oname == 'pk4-diff' and truth_data.get('pk4-sample-heads') is not None:
            df = truth_data['pk4-sample-heads']
            if kper == 1:
                time_idx = 1
            elif kper == 2:
                time_idx = 1 + kstp
            elif kper == 3:
                time_idx = 54
            elif kper == 4:
                time_idx = 54 + kstp
            else:
                return None

            if 'time' in df.columns and 'pk4-diff' in df.columns:
                row = df[df['time'] == time_idx]
                if len(row) > 0:
                    return row['pk4-diff'].values[0]
            return None

        # pk4-spring-diff - use pk4-spr-diff column (backward compatibility)
        elif oname == 'pk4-spring-diff' and truth_data.get('pk4-sample-heads') is not None:
            df = truth_data['pk4-sample-heads']
            if kper == 1:
                time_idx = 1
            elif kper == 2:
                time_idx = 1 + kstp
            elif kper == 3:
                time_idx = 54
            elif kper == 4:
                time_idx = 54 + kstp
            else:
                return None

            if 'time' in df.columns and 'pk4-spr-diff' in df.columns:
                row = df[df['time'] == time_idx]
                if len(row) > 0:
                    return row['pk4-spr-diff'].values[0]
            return None

        # pk4-aw-diff - use pk4-aw-diff column (backward compatibility)
        elif oname == 'pk4-aw-diff' and truth_data.get('pk4-sample-heads') is not None:
            df = truth_data['pk4-sample-heads']
            if kper == 1:
                time_idx = 1
            elif kper == 2:
                time_idx = 1 + kstp
            elif kper == 3:
                time_idx = 54
            elif kper == 4:
                time_idx = 54 + kstp
            else:
                return None

            if 'time' in df.columns and 'pk4-aw-diff' in df.columns:
                row = df[df['time'] == time_idx]
                if len(row) > 0:
                    return row['pk4-aw-diff'].values[0]
            return None

        # pk4-heads
        elif oname == 'pk4-heads' and truth_data.get('pk4-heads') is not None:
            df = truth_data['pk4-heads']
            # Map kper/kstp to time index
            if kper == 1:
                time_idx = 1
            elif kper == 2:
                time_idx = 1 + kstp
            elif kper == 3:
                time_idx = 54
            elif kper == 4:
                time_idx = 54 + kstp
            else:
                return None

            # Find row with this time
            if 'time' in df.columns and usecol in df.columns:
                row = df[df['time'] == time_idx]
                if len(row) > 0:
                    return row[usecol].values[0]
            return None

        # budget-sw - use sw column
        elif oname == 'budget-sw' and truth_data.get('budget') is not None:
            df = truth_data['budget']
            if 'kper' in df.columns and 'sw' in df.columns:
                row = df[df['kper'] == kper]
                if len(row) > 0:
                    return row['sw'].values[0]
            return None

        # budget-rch - use rch or recharge column
        elif oname == 'budget-rch' and truth_data.get('budget') is not None:
            df = truth_data['budget']
            if 'kper' in df.columns:
                row = df[df['kper'] == kper]
                if len(row) > 0:
                    # Try 'rch' first, then 'recharge'
                    if 'rch' in df.columns:
                        return row['rch'].values[0]
                    elif 'recharge' in df.columns:
                        return row['recharge'].values[0]
            return None

        # budget (original handler for other usecol values)
        elif oname == 'budget' and truth_data.get('budget') is not None:
            df = truth_data['budget']
            if 'kper' in df.columns:
                row = df[df['kper'] == kper]
                if len(row) > 0:
                    # Try usecol first, then check for common aliases
                    if usecol in df.columns:
                        return row[usecol].values[0]
                    elif usecol == 'rch' and 'recharge' in df.columns:
                        return row['recharge'].values[0]
                    elif usecol == 'recharge' and 'rch' in df.columns:
                        return row['rch'].values[0]
            return None

        # recession - match by kper and usecol
        elif oname == 'recession' and truth_data.get('recession') is not None:
            df = truth_data['recession']
            if 'kper' in df.columns and usecol in df.columns:
                row = df[df['kper'] == kper]
                if len(row) > 0:
                    return row[usecol].values[0]
            return None

        # aw-q - use AWq column from GHB_AW_fluxes
        elif oname == 'aw-q' and truth_data.get('awq') is not None:
            df = truth_data['awq']
            i = obs_row.get('i')
            j = obs_row.get('j')
            if i is not None and j is not None:
                i = int(i)
                j = int(j)
                # Match by kper, kstp, i, j
                mask = (df['kper'] == kper) & (df['kstp'] == kstp) & (df['i'] == i) & (df['j'] == j)
                row = df[mask]
                if len(row) > 0 and 'AWq' in df.columns:
                    return row['AWq'].values[0]
            return None

        # arr-awq - match by kper, kstp, i, j (original handler)
        elif oname == 'arr-awq' and truth_data.get('awq') is not None:
            df = truth_data['awq']
            i = obs_row.get('i')
            j = obs_row.get('j')
            if i is not None and j is not None:
                i = int(i)
                j = int(j)
                # Match by kper, kstp, i, j
                mask = (df['kper'] == kper) & (df['kstp'] == kstp) & (df['i'] == i) & (df['j'] == j)
                row = df[mask]
                if len(row) > 0 and 'AWq' in df.columns:
                    return row['AWq'].values[0]
            return None

        # spring-q - use SPq column from GHB_SP_fluxes
        elif oname == 'spring-q' and truth_data.get('spq') is not None:
            df = truth_data['spq']
            i = obs_row.get('i')
            j = obs_row.get('j')
            if i is not None and j is not None:
                i = int(i)
                j = int(j)
                # Match by kper, kstp, i, j
                mask = (df['kper'] == kper) & (df['kstp'] == kstp) & (df['i'] == i) & (df['j'] == j)
                row = df[mask]
                if len(row) > 0 and 'SPq' in df.columns:
                    return row['SPq'].values[0]
            return None

        # arr-spq - match by kper, kstp, i, j (original handler)
        elif oname == 'arr-spq' and truth_data.get('spq') is not None:
            df = truth_data['spq']
            i = obs_row.get('i')
            j = obs_row.get('j')
            if i is not None and j is not None:
                i = int(i)
                j = int(j)
                # Match by kper, kstp, i, j
                mask = (df['kper'] == kper) & (df['kstp'] == kstp) & (df['i'] == i) & (df['j'] == j)
                row = df[mask]
                if len(row) > 0 and 'SPq' in df.columns:
                    return row['SPq'].values[0]
            return None

        # top-heads - excluded from plots, but handler provided for completeness
        elif oname == 'top-heads' and truth_data.get('top-heads') is not None:
            df = truth_data['top-heads']
            i = obs_row.get('i')
            j = obs_row.get('j')
            if i is not None and j is not None:
                i = int(i)
                j = int(j)
                # Match by kper, kstp, i, j
                # Assuming the .dat file has columns: kper, kstp, i, j, head
                if all(col in df.columns for col in ['kper', 'kstp', 'i', 'j']):
                    mask = (df['kper'] == kper) & (df['kstp'] == kstp) & (df['i'] == i) & (df['j'] == j)
                    row = df[mask]
                    if len(row) > 0:
                        # Try common column names for head values
                        for head_col in ['head', 'heads', 'value']:
                            if head_col in df.columns:
                                return row[head_col].values[0]
            return None

        else:
            return None

    except Exception as e:
        return None


def is_constraint_obs(obgnme):
    """Check if observation is a constraint (less_than or greater_than)."""
    obgnme_str = str(obgnme).lower()
    return obgnme_str.startswith('less_') or obgnme_str.startswith('greater_')


def create_1to1_plots(run_name, model_name, post_iter=19, filter_file=None, suffix=None):
    """
    Create 1:1 plots for each observation group and usecol combination.
    """
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"1:1 PLOTS BY OBSERVATION GROUP: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    master_dir = os.path.join('models', run_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', run_name, 'pp_1to1_plots')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data.copy()

    print(f"\nTotal observations: {len(obs_data)}")

    # Load truth data
    print("\nLoading truth data...")
    truth_data = load_truth_data(run_name)

    # Filter observations:
    # 1. Non-zero weight
    # 2. Not a constraint (less_than/greater_than)
    # 3. Not top-heads (excluded from plots per user requirement)
    obs_data['is_constraint'] = obs_data['obgnme'].apply(is_constraint_obs)
    obs_data['is_top_heads'] = obs_data['oname'].apply(lambda x: str(x).lower() == 'top-heads')
    filtered_obs = obs_data[(obs_data['weight'] > 0) &
                            (~obs_data['is_constraint']) &
                            (~obs_data['is_top_heads'])].copy()

    print(f"\nFiltered observations (non-zero weight, non-constraint, excluding top-heads): {len(filtered_obs)}")

    # Group by oname and usecol
    if 'usecol' not in filtered_obs.columns:
        filtered_obs['usecol'] = ''

    # Get unique oname/usecol combinations
    groups = filtered_obs.groupby(['oname', 'usecol']).size().reset_index(name='count')
    groups = groups[groups['count'] > 0]

    print(f"\nObservation groups to plot:")
    for _, row in groups.iterrows():
        print(f"  {row['oname']} / {row['usecol']}: {row['count']} observations")

    # Load filtered realizations
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations")

    # Load posterior observation ensemble
    print(f"\nLoading posterior observation ensemble (iter {post_iter})...")
    post_file = os.path.join(master_dir, f'{model_name}.{post_iter}.obs.csv')
    post_oe = pd.read_csv(post_file, index_col=0, low_memory=False)

    # Apply filter
    if filtered_ids:
        post_oe = post_oe.loc[[r for r in filtered_ids if r in post_oe.index]]
        print(f"  Filtered to {len(post_oe)} realizations")

    # Prepare data for each group
    plot_data = []

    for _, group_row in groups.iterrows():
        oname = group_row['oname']
        usecol = group_row['usecol']

        # Get observations for this group
        if usecol:
            group_obs = filtered_obs[(filtered_obs['oname'] == oname) &
                                     (filtered_obs['usecol'] == usecol)]
        else:
            group_obs = filtered_obs[filtered_obs['oname'] == oname]

        if len(group_obs) == 0:
            continue

        # Get truth and simulated values
        truth_vals = []
        sim_vals = []
        sim_min = []  # Actual min for plotting
        sim_max = []  # Actual max for plotting
        sim_min_clipped = []  # IQR-clipped min for axis limits
        sim_max_clipped = []  # IQR-clipped max for axis limits
        sim_p25 = []
        sim_p75 = []
        obs_names = []

        for obs_name, obs_row in group_obs.iterrows():
            # Get truth value
            truth_val = get_truth_value(obs_row, truth_data)

            # Get simulated values (percentiles of posterior ensemble)
            if obs_name in post_oe.columns:
                ensemble_vals = post_oe[obs_name].values
                sim_val = np.median(ensemble_vals)
                actual_min = np.min(ensemble_vals)
                actual_max = np.max(ensemble_vals)
                p25 = np.percentile(ensemble_vals, 25)
                p75 = np.percentile(ensemble_vals, 75)
                IQR = p75 - p25
                # Clip to 1.5*IQR for axis limits
                clipped_min = max(actual_min, p25 - 1.5 * IQR)
                clipped_max = min(actual_max, p75 + 1.5 * IQR)

                
            else:
                sim_val = None

            if truth_val is not None and sim_val is not None:
                truth_vals.append(truth_val)
                sim_vals.append(sim_val)
                sim_min.append(actual_min)
                sim_max.append(actual_max)
                sim_min_clipped.append(clipped_min)
                sim_max_clipped.append(clipped_max)
                sim_p25.append(p25)
                sim_p75.append(p75)
                obs_names.append(obs_name)

        if len(truth_vals) > 0:
            plot_data.append({
                'oname': oname,
                'usecol': usecol,
                'truth': np.array(truth_vals),
                'simulated': np.array(sim_vals),
                'sim_min': np.array(sim_min),
                'sim_max': np.array(sim_max),
                'sim_min_clipped': np.array(sim_min_clipped),
                'sim_max_clipped': np.array(sim_max_clipped),
                'sim_p25': np.array(sim_p25),
                'sim_p75': np.array(sim_p75),
                'obs_names': obs_names,
                'label': f"{oname}:{usecol}" if usecol else oname
            })

    print(f"\nGroups with data: {len(plot_data)}")

    if len(plot_data) == 0:
        print("No data to plot!")
        return

    # Create figure with grid layout (3 columns)
    n_plots = len(plot_data)
    n_cols = 3
    n_rows = int(np.ceil(n_plots / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))

    # Flatten axes for easy iteration
    if n_rows == 1 and n_cols == 1:
        axes = np.array([[axes]])
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)

    axes_flat = axes.flatten()

    # Plot each group
    for idx, data in enumerate(plot_data):
        ax = axes_flat[idx]

        truth = data['truth']
        sim = data['simulated']
        sim_min = data['sim_min']
        sim_p25 = data['sim_p25']
        sim_p75 = data['sim_p75']
        sim_max = data['sim_max']
        sim_min_clipped = data['sim_min_clipped']
        sim_max_clipped = data['sim_max_clipped']
        label = data['label']

        # Set axis limits based on IQR-clipped values for both axes
        sim_all = np.concatenate([sim, sim_min_clipped, sim_max_clipped])
        sim_min_val, sim_max_val = np.min(sim_all), np.max(sim_all)
        margin = (sim_max_val - sim_min_val) * 0.1 if sim_max_val > sim_min_val else 0.1
        axis_range = [sim_min_val - margin, sim_max_val + margin]
        
        # Plot uncertainty ranges for all points
        # Min to max range (blue, thick)
        ax.vlines(truth, sim_min, sim_max, colors="#000000", linewidth=0.5, alpha=1, zorder=1)
        # 25th to 75th percentile (black, thicker on top) - use loop to ensure all points
        for i in range(len(truth)):
            ax.plot([truth[i], truth[i]], [sim_p25[i], sim_p75[i]],
                   color='black', linewidth=8, alpha=0.5, zorder=2)

        # Plot median as small horizontal dark red line
        # Calculate line half-width based on data range
        x_range = axis_range[1] - axis_range[0]
        line_half_width = x_range * 0.015
        for i in range(len(truth)):
            ax.plot([truth[i] - line_half_width, truth[i] + line_half_width],
                   [sim[i], sim[i]], color='darkred', linewidth=1, zorder=3)

        # Manual adjustment for arr-awq
        if data['oname'] == 'arr-awq':
            axis_range[1] = -400  # Set max to -400

        # Plot 1:1 line
        ax.plot(axis_range, axis_range, 'k--', linewidth=1, alpha=0.5, label='1:1')

        ax.set_xlabel('Truth', fontsize=8)
        ax.set_ylabel('Simulated (median)', fontsize=8)
        ax.set_title(label, fontsize=9, fontweight='bold')
        ax.set_xlim(axis_range)
        ax.set_ylim(axis_range)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)

    # Hide empty subplots
    for idx in range(len(plot_data), len(axes_flat)):
        axes_flat[idx].set_visible(False)

    plt.tight_layout()

    # Save figure
    output_path = os.path.join(output_dir, f'{run_name}_1to1_plots_by_obsgroup{file_suffix}.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"\nSaved plot to: {output_path}")

    # Save summary CSV with individual observations
    summary_rows = []
    for data in plot_data:
        truth = data['truth']
        sim = data['simulated']
        obs_names = data['obs_names']

        for i, obs_name in enumerate(obs_names):
            truth_val = truth[i]
            sim_val = sim[i]
            residual = sim_val - truth_val

            # Calculate normalized RMSE (as percentage of truth range for this group)
            truth_range = np.max(truth) - np.min(truth) if len(truth) > 1 else 1
            nrmse = (abs(residual) / truth_range * 100) if truth_range > 0 else 0

            summary_rows.append({
                'oname': data['oname'],
                'usecol': data['usecol'],
                'obs_name': obs_name,
                'truth': truth_val,
                'simulated': sim_val,
                'residual': residual,
                'nrmse_pct': nrmse,
            })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = os.path.join(output_dir, f'{run_name}_1to1_summary{file_suffix}.csv')
    summary_df.to_csv(summary_path, index=False)
    print(f"Saved summary to: {summary_path}")

    # Save aggregated statistics by oname and kper
    kper_stats_rows = []

    # Get all unique onames
    unique_onames = filtered_obs['oname'].unique()

    for oname in unique_onames:
        # Get observations for this oname
        oname_obs = filtered_obs[filtered_obs['oname'] == oname].copy()

        # Convert kper to numeric
        oname_obs['kper'] = pd.to_numeric(oname_obs['kper'], errors='coerce')

        # Get unique kpers
        unique_kpers = sorted(oname_obs['kper'].dropna().unique())

        for kper in unique_kpers:
            kper = int(kper)
            # Get observations for this oname/kper
            kper_obs = oname_obs[oname_obs['kper'] == kper]

            if len(kper_obs) == 0:
                continue

            # Collect truth and simulated values
            truth_vals = []
            sim_vals = []

            for obs_name, obs_row in kper_obs.iterrows():
                # Get truth value
                truth_val = get_truth_value(obs_row, truth_data)

                # Get simulated median
                if obs_name in post_oe.columns:
                    sim_val = np.median(post_oe[obs_name].values)
                else:
                    sim_val = None

                if truth_val is not None and sim_val is not None:
                    truth_vals.append(truth_val)
                    sim_vals.append(sim_val)

            if len(truth_vals) > 0:
                truth_arr = np.array(truth_vals)
                sim_arr = np.array(sim_vals)
                residuals = sim_arr - truth_arr

                kper_stats_rows.append({
                    'oname': oname,
                    'kper': kper,
                    'n_obs': len(truth_vals),
                    'simulated_median': np.median(sim_arr),
                    'observed_median': np.median(truth_arr),
                    'rmse': np.sqrt(np.mean(residuals**2)),
                })

    if kper_stats_rows:
        kper_stats_df = pd.DataFrame(kper_stats_rows)
        kper_stats_path = os.path.join(output_dir, f'{run_name}_1to1_kper_stats{file_suffix}.csv')
        kper_stats_df.to_csv(kper_stats_path, index=False)
        print(f"Saved kper stats to: {kper_stats_path}")

    # Save group statistics
    stats_rows = []
    for data in plot_data:
        truth = data['truth']
        sim = data['simulated']
        residuals = sim - truth

        rmse = np.sqrt(np.mean(residuals**2))
        truth_range = np.max(truth) - np.min(truth) if len(truth) > 1 else 1
        nrmse = (rmse / truth_range * 100) if truth_range > 0 else 0

        stats_rows.append({
            'oname': data['oname'],
            'usecol': data['usecol'],
            'n_obs': len(truth),
            'rmse': rmse,
            'nrmse_pct': nrmse,
            'mae': np.mean(np.abs(residuals)),
            'r2': 1 - np.sum(residuals**2) / np.sum((truth - np.mean(truth))**2) if len(truth) > 1 else 0,
            'min_truth': np.min(truth),
            'max_truth': np.max(truth),
            'min_sim': np.min(sim),
            'max_sim': np.max(sim),
        })

    stats_df = pd.DataFrame(stats_rows)
    stats_path = os.path.join(output_dir, f'{run_name}_1to1_stats{file_suffix}.csv')
    stats_df.to_csv(stats_path, index=False)
    print(f"Saved stats to: {stats_path}")

    print(f"\nAll outputs saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_1to1_plots(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        suffix=args.suffix
    )
