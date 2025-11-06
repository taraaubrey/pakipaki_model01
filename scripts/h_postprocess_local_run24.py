"""
Comprehensive post-processing script for local_run24 master_ies results.

This script generates:
1. Phi progression boxplots through iterations
2. Prior/posterior observation ensemble comparisons (grouped by obsgnme)
   - Special handling for array observations (i,j coordinates) grouped by kper/kstp
3. Prior/posterior parameter ensemble comparisons (grouped by pargnme)

Usage:
    python scripts/h_postprocess_local_run24.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path
import pyemu

# Configuration
MODEL_NAME = 'local_run24'
MASTER_IES_DIR = f'models/{MODEL_NAME}/pest/master_ies'
OUTPUT_DIR = f'models/{MODEL_NAME}/figures'

# Create output directory
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Helper functions to parse observation names
def extract_kper(name):
    """Extract kper from observation name"""
    match = re.search(r'kper(\d+)', name)
    return int(match.group(1)) if match else None

def extract_kstp(name):
    """Extract kstp from observation name"""
    match = re.search(r'kstp(\d+)', name)
    return int(match.group(1)) if match else None

def extract_i(name):
    """Extract i index from observation name"""
    match = re.search(r'_i:(\d+)', name)
    return int(match.group(1)) if match else None

def extract_j(name):
    """Extract j index from observation name"""
    match = re.search(r'_j:(\d+)', name)
    return int(match.group(1)) if match else None

def has_ij_coords(name):
    """Check if observation name has i,j coordinates"""
    return ('_i:' in name) and ('_j:' in name)


def plot_phi_progression():
    """
    Plot phi progression through iterations as boxplots.
    """
    print("\n1. Plotting phi progression...")

    phi_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.phi.actual.csv')
    if not os.path.exists(phi_file):
        print(f"  WARNING: Phi file not found: {phi_file}")
        return

    # Read phi data
    phi_df = pd.read_csv(phi_file, index_col=0)

    # Get iterations and ensemble columns
    iterations = phi_df.index.tolist()
    # Columns after 'max' are individual realizations
    realization_cols = phi_df.columns[6:]  # Skip iteration, total_runs, mean, std, min, max

    # Prepare data for boxplot
    phi_data = []
    phi_labels = []
    for iteration in iterations:
        phi_values = phi_df.loc[iteration, realization_cols].values.astype(float)
        phi_data.append(phi_values)
        phi_labels.append(f'Iter {iteration}')

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 6))

    # Create boxplot
    bp = ax.boxplot(phi_data, labels=phi_labels, patch_artist=True)

    # Customize boxplot colors
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
        patch.set_alpha(0.7)

    ax.set_yscale('log')
    ax.set_ylabel('Phi (log scale)', fontsize=12)
    ax.set_xlabel('Iteration', fontsize=12)
    ax.set_title(f'Phi Progression Through IES Iterations - {MODEL_NAME}', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')

    plt.xticks(rotation=45)
    plt.tight_layout()

    output_file = os.path.join(OUTPUT_DIR, 'phi_progression_boxplot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")


def plot_observation_ensembles():
    """
    Plot prior and posterior observation ensembles grouped by obsgnme.
    For array observations (with i,j), create spatial plots grouped by kper/kstp.
    """
    print("\n2. Plotting observation ensembles...")

    # Load PST to get observation groups
    pst_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.pst')
    pst = pyemu.Pst(pst_file)

    # Find prior and posterior files
    prior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.0.obs.csv')

    # Find last iteration
    obs_files = sorted(Path(MASTER_IES_DIR).glob(f'{MODEL_NAME}.*.obs.csv'))
    iter_numbers = [int(f.stem.split('.')[-2]) for f in obs_files if f.stem.split('.')[-2].isdigit()]
    last_iter = max(iter_numbers)
    posterior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.{last_iter}.obs.csv')

    print(f"  Prior: iteration 0")
    print(f"  Posterior: iteration {last_iter}")

    if not os.path.exists(prior_file) or not os.path.exists(posterior_file):
        print(f"  ERROR: Observation files not found")
        return

    # Read observation data (these files are large, so we'll read in chunks)
    print(f"  Loading observation ensembles...")
    prior_df = pd.read_csv(prior_file, index_col=0)
    posterior_df = pd.read_csv(posterior_file, index_col=0)

    # Get observation groups
    obs_groups = pst.observation_data.groupby('obgnme').groups

    print(f"  Found {len(obs_groups)} observation groups")

    # Create output subdirectory
    obs_output_dir = os.path.join(OUTPUT_DIR, 'observations')
    os.makedirs(obs_output_dir, exist_ok=True)

    # Process each observation group
    for obgnme, obs_indices in obs_groups.items():
        obs_names = pst.observation_data.loc[obs_indices].index.tolist()

        # Check if these are array observations
        if any(has_ij_coords(name) for name in obs_names):
            print(f"\n  Processing array observations for group: {obgnme}")
            plot_array_observations(obgnme, obs_names, prior_df, posterior_df, obs_output_dir)
        else:
            print(f"\n  Processing time series observations for group: {obgnme}")
            plot_timeseries_obs_histograms(obgnme, obs_names, prior_df, posterior_df, obs_output_dir)


def plot_array_observations(obgnme, obs_names, prior_df, posterior_df, output_dir):
    """
    Plot array observations grouped by kper and kstp.
    Creates 3-column plots: mean, min, max
    """
    # Group observations by kper and kstp
    obs_data = []
    for name in obs_names:
        if name not in prior_df.columns:
            continue
        obs_data.append({
            'name': name,
            'kper': extract_kper(name),
            'kstp': extract_kstp(name),
            'i': extract_i(name),
            'j': extract_j(name)
        })

    obs_df = pd.DataFrame(obs_data)

    if obs_df.empty or obs_df['kper'].isna().all():
        print(f"    No valid array observations with kper/kstp found")
        return

    # Group by kper and kstp
    grouped = obs_df.groupby(['kper', 'kstp'])

    print(f"    Found {len(grouped)} unique kper/kstp combinations")

    # Get grid extent
    max_i = obs_df['i'].max()
    max_j = obs_df['j'].max()

    # Process each kper/kstp combination
    for (kper, kstp), group in grouped:
        if pd.isna(kper) or pd.isna(kstp):
            continue

        kper = int(kper)
        kstp = int(kstp)

        # Initialize arrays
        prior_mean = np.full((max_i + 1, max_j + 1), np.nan)
        prior_min = np.full((max_i + 1, max_j + 1), np.nan)
        prior_max = np.full((max_i + 1, max_j + 1), np.nan)

        post_mean = np.full((max_i + 1, max_j + 1), np.nan)
        post_min = np.full((max_i + 1, max_j + 1), np.nan)
        post_max = np.full((max_i + 1, max_j + 1), np.nan)

        # Fill arrays
        for _, row in group.iterrows():
            name = row['name']
            i = int(row['i'])
            j = int(row['j'])

            if name in prior_df.columns:
                prior_vals = prior_df[name].values
                prior_mean[i, j] = np.mean(prior_vals)
                prior_min[i, j] = np.min(prior_vals)
                prior_max[i, j] = np.max(prior_vals)

            if name in posterior_df.columns:
                post_vals = posterior_df[name].values
                post_mean[i, j] = np.mean(post_vals)
                post_min[i, j] = np.min(post_vals)
                post_max[i, j] = np.max(post_vals)

        # Create figure with 2 rows (prior, posterior) x 3 columns (mean, min, max)
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        fig.suptitle(f'{obgnme} - kper={kper}, kstp={kstp}', fontsize=14, fontweight='bold')

        # Determine common colorbar limits
        all_data = np.concatenate([
            prior_mean[~np.isnan(prior_mean)],
            post_mean[~np.isnan(post_mean)]
        ])
        if len(all_data) > 0:
            vmin, vmax = np.percentile(all_data, [2, 98])
        else:
            vmin, vmax = 0, 1

        # Prior row
        im = axes[0, 0].imshow(prior_mean, cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
        axes[0, 0].set_title('Prior - Mean')
        axes[0, 0].set_ylabel('i index')
        plt.colorbar(im, ax=axes[0, 0])

        im = axes[0, 1].imshow(prior_min, cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
        axes[0, 1].set_title('Prior - Min')
        plt.colorbar(im, ax=axes[0, 1])

        im = axes[0, 2].imshow(prior_max, cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
        axes[0, 2].set_title('Prior - Max')
        plt.colorbar(im, ax=axes[0, 2])

        # Posterior row
        im = axes[1, 0].imshow(post_mean, cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
        axes[1, 0].set_title('Posterior - Mean')
        axes[1, 0].set_ylabel('i index')
        axes[1, 0].set_xlabel('j index')
        plt.colorbar(im, ax=axes[1, 0])

        im = axes[1, 1].imshow(post_min, cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
        axes[1, 1].set_title('Posterior - Min')
        axes[1, 1].set_xlabel('j index')
        plt.colorbar(im, ax=axes[1, 1])

        im = axes[1, 2].imshow(post_max, cmap='viridis', origin='upper', vmin=vmin, vmax=vmax)
        axes[1, 2].set_title('Posterior - Max')
        axes[1, 2].set_xlabel('j index')
        plt.colorbar(im, ax=axes[1, 2])

        plt.tight_layout()

        output_file = os.path.join(output_dir, f'{obgnme}_kper{kper}_kstp{kstp}.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {output_file}")


def plot_timeseries_obs_histograms(obgnme, obs_names, prior_df, posterior_df, output_dir):
    """
    Plot time series observations as histograms comparing prior and posterior.
    """
    # Filter to observations that exist in both prior and posterior
    valid_obs = [name for name in obs_names if name in prior_df.columns and name in posterior_df.columns]

    if len(valid_obs) == 0:
        print(f"    No valid observations found for group {obgnme}")
        return

    # Limit number of observations to plot (for large groups)
    max_obs_per_plot = 20
    if len(valid_obs) > max_obs_per_plot:
        print(f"    Warning: {len(valid_obs)} observations in group, plotting first {max_obs_per_plot}")
        valid_obs = valid_obs[:max_obs_per_plot]

    # Create figure
    n_obs = len(valid_obs)
    ncols = min(4, n_obs)
    nrows = int(np.ceil(n_obs / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows))
    fig.suptitle(f'Observation Ensembles - {obgnme}', fontsize=14, fontweight='bold')

    if n_obs == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, obs_name in enumerate(valid_obs):
        ax = axes[idx]

        # Get prior and posterior values
        prior_vals = prior_df[obs_name].values
        post_vals = posterior_df[obs_name].values

        # Plot histograms
        ax.hist(prior_vals, bins=30, alpha=0.5, label='Prior', color='blue')
        ax.hist(post_vals, bins=30, alpha=0.5, label='Posterior', color='red')

        ax.set_title(obs_name[:40], fontsize=8)  # Truncate long names
        ax.set_xlabel('Value')
        ax.set_ylabel('Frequency')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    # Hide unused subplots
    for idx in range(n_obs, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()

    output_file = os.path.join(output_dir, f'{obgnme}_histograms.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {output_file}")


def plot_parameter_ensembles():
    """
    Plot prior and posterior parameter ensembles grouped by pargnme.
    """
    print("\n3. Plotting parameter ensembles...")

    # Load PST to get parameter groups
    pst_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.pst')
    pst = pyemu.Pst(pst_file)

    # Find prior and posterior files
    prior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.0.par.csv')

    # Find last iteration
    par_files = sorted(Path(MASTER_IES_DIR).glob(f'{MODEL_NAME}.*.par.csv'))
    iter_numbers = [int(f.stem.split('.')[-2]) for f in par_files if f.stem.split('.')[-2].isdigit()]
    last_iter = max(iter_numbers)
    posterior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.{last_iter}.par.csv')

    print(f"  Prior: iteration 0")
    print(f"  Posterior: iteration {last_iter}")

    if not os.path.exists(prior_file) or not os.path.exists(posterior_file):
        print(f"  ERROR: Parameter files not found")
        return

    # Read parameter data
    print(f"  Loading parameter ensembles...")
    prior_df = pd.read_csv(prior_file, index_col=0)
    posterior_df = pd.read_csv(posterior_file, index_col=0)

    # Get parameter groups
    par_groups = pst.parameter_data.groupby('pargp').groups

    print(f"  Found {len(par_groups)} parameter groups")

    # Create output subdirectory
    par_output_dir = os.path.join(OUTPUT_DIR, 'parameters')
    os.makedirs(par_output_dir, exist_ok=True)

    # Process each parameter group
    for pargnme, par_indices in par_groups.items():
        par_names = pst.parameter_data.loc[par_indices].index.tolist()

        # Filter to parameters that exist in both prior and posterior
        valid_pars = [name for name in par_names if name in prior_df.columns and name in posterior_df.columns]

        if len(valid_pars) == 0:
            continue

        print(f"\n  Processing parameter group: {pargnme} ({len(valid_pars)} parameters)")

        # Limit number of parameters to plot
        max_pars_per_plot = 20
        if len(valid_pars) > max_pars_per_plot:
            print(f"    Warning: {len(valid_pars)} parameters in group, plotting first {max_pars_per_plot}")
            valid_pars = valid_pars[:max_pars_per_plot]

        # Create figure
        n_pars = len(valid_pars)
        ncols = min(4, n_pars)
        nrows = int(np.ceil(n_pars / ncols))

        fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows))
        fig.suptitle(f'Parameter Ensembles - {pargnme}', fontsize=14, fontweight='bold')

        if n_pars == 1:
            axes = np.array([axes])
        axes = axes.flatten()

        for idx, par_name in enumerate(valid_pars):
            ax = axes[idx]

            # Get prior and posterior values
            prior_vals = prior_df[par_name].values
            post_vals = posterior_df[par_name].values

            # Plot histograms
            ax.hist(prior_vals, bins=30, alpha=0.5, label='Prior', color='blue')
            ax.hist(post_vals, bins=30, alpha=0.5, label='Posterior', color='red')

            ax.set_title(par_name[:40], fontsize=8)  # Truncate long names
            ax.set_xlabel('Value')
            ax.set_ylabel('Frequency')
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        # Hide unused subplots
        for idx in range(n_pars, len(axes)):
            axes[idx].set_visible(False)

        plt.tight_layout()

        output_file = os.path.join(par_output_dir, f'{pargnme}_histograms.png')
        plt.savefig(output_file, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {output_file}")


def extract_time(name):
    """Extract time index from observation name"""
    match = re.search(r'time:(\d+)', name)
    return int(match.group(1)) if match else None

def extract_usecol(name):
    """Extract usecol from observation name"""
    match = re.search(r'usecol:([^_]+)', name)
    return match.group(1) if match else None


def plot_timeseries_observations():
    """
    Plot time series observations comparing prior and posterior with histograms and truth data.

    Data structure: rows = realizations, columns = observations
    """
    print("\n4. Plotting time series observations (heads)...")

    # Load PST to get headers
    pst_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.pst')
    pst = pyemu.Pst(pst_file)
    headers = pst.observation_data.index.tolist()

    # Find prior and posterior files
    prior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.0.obs.csv')
    par_files = sorted(Path(MASTER_IES_DIR).glob(f'{MODEL_NAME}.*.obs.csv'))
    iter_numbers = [int(f.stem.split('.')[-2]) for f in par_files if f.stem.split('.')[-2].isdigit()]
    last_iter = max(iter_numbers)
    posterior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.{last_iter}.obs.csv')

    # Look for truth data
    truth_file = f'models/{MODEL_NAME}/{MODEL_NAME}/truth/output.sample_heads.truth.csv'
    if os.path.exists(truth_file):
        print(f"  Loading truth data from: {truth_file}")
        truth_df = pd.read_csv(truth_file, index_col=0)
    else:
        print(f"  No truth data found at: {truth_file}")
        truth_df = None

    # Look for time series head observations
    col_startswith = 'ts-heads'
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("  No time series head observations found, skipping")
        return

    print(f"  Loading {len(desired_cols)} time series observations...")

    # Read data - rows are realizations, columns are observations
    pr_oe_ts = pd.read_csv(prior_file, usecols=desired_cols, index_col=0)
    pt_oe_ts = pd.read_csv(posterior_file, usecols=desired_cols, index_col=0)

    # Extract metadata from column names
    col_metadata = []
    for col in desired_cols:
        time_val = extract_time(col)
        usecol_val = extract_usecol(col)
        col_metadata.append({'col': col, 'time': time_val, 'usecol': usecol_val})

    metadata_df = pd.DataFrame(col_metadata)

    # Get unique locations
    unique_usecols = metadata_df['usecol'].unique()
    n = len(unique_usecols)

    print(f"  Creating plots for {n} time series locations...")

    # Create figure with gridspec (time series + histogram for each location)
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(14, n * 3))
    gs = gridspec.GridSpec(n, 2, width_ratios=[4, 1], figure=fig)

    for idx, usecol in enumerate(unique_usecols):
        # Create subplot axes
        ax_ts = fig.add_subplot(gs[idx, 0])
        ax_hist = fig.add_subplot(gs[idx, 1], sharey=ax_ts)

        # Get columns for this location
        usecol_cols = metadata_df[metadata_df['usecol'] == usecol].sort_values('time')
        time_values = usecol_cols['time'].values
        obs_cols = usecol_cols['col'].values

        # Get data for this location from both prior and posterior
        prior_data = pr_oe_ts[obs_cols].replace(-999, np.nan)  # shape: (n_realizations, n_times)
        post_data = pt_oe_ts[obs_cols].replace(-999, np.nan)   # shape: (n_realizations, n_times)

        # Calculate percentiles for posterior (1st to 99th) across all times and realizations
        all_post_values = post_data.values.flatten()
        all_post_values = all_post_values[~np.isnan(all_post_values)]

        if len(all_post_values) > 0:
            y_min = np.percentile(all_post_values, 1)
            y_max = np.percentile(all_post_values, 99)
        else:
            y_min, y_max = None, None

        # Plot prior ensemble (grey) - each row is a realization
        for real_idx in prior_data.index:
            y_vals = prior_data.loc[real_idx].values
            # Only plot if within percentile range
            if y_min is not None and y_max is not None:
                # Check if any values are in range
                if np.any((y_vals >= y_min) & (y_vals <= y_max) & ~np.isnan(y_vals)):
                    ax_ts.plot(time_values, y_vals, alpha=0.2, color='grey', linewidth=0.5)
            else:
                ax_ts.plot(time_values, y_vals, alpha=0.2, color='grey', linewidth=0.5)

        # Plot posterior ensemble (blue) - each row is a realization
        for real_idx in post_data.index:
            y_vals = post_data.loc[real_idx].values
            # Only plot if within percentile range
            if y_min is not None and y_max is not None:
                # Check if any values are in range
                if np.any((y_vals >= y_min) & (y_vals <= y_max) & ~np.isnan(y_vals)):
                    ax_ts.plot(time_values, y_vals, alpha=0.3, color='blue', linewidth=0.8)
            else:
                ax_ts.plot(time_values, y_vals, alpha=0.3, color='blue', linewidth=0.8)

        # Plot truth data if available
        if truth_df is not None and usecol in truth_df.columns:
            truth_vals = truth_df[usecol].values
            ax_ts.plot(range(len(truth_vals)), truth_vals, 'r-', linewidth=2, label='Truth', zorder=10)

        # Set y-limits based on posterior 1-99 percentiles
        if y_min is not None and y_max is not None:
            y_range = y_max - y_min
            ax_ts.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

        ax_ts.set_title(f'{usecol}', loc='left', fontsize=10)
        ax_ts.set_ylabel('Head (m)')
        ax_ts.grid(alpha=0.3)

        if idx == 0:
            legend_items = ['Prior', 'Posterior']
            if truth_df is not None:
                legend_items.append('Truth')
            ax_ts.legend(legend_items, fontsize=8)

        if idx == n - 1:
            ax_ts.set_xlabel('Time')

        # RIGHT PLOT: Histogram
        # Collect all values for histograms (all times, all realizations)
        prior_all_vals = prior_data.values.flatten()
        prior_all_vals = prior_all_vals[~np.isnan(prior_all_vals)]

        post_all_vals = post_data.values.flatten()
        post_all_vals = post_all_vals[~np.isnan(post_all_vals)]

        # Filter to 1-99 percentile range
        if y_min is not None and y_max is not None:
            prior_filtered = prior_all_vals[(prior_all_vals >= y_min) & (prior_all_vals <= y_max)]
            post_filtered = post_all_vals[(post_all_vals >= y_min) & (post_all_vals <= y_max)]
        else:
            prior_filtered = prior_all_vals
            post_filtered = post_all_vals

        # Plot histograms horizontally
        if len(prior_filtered) > 0:
            ax_hist.hist(prior_filtered, bins=30, alpha=0.5, color='grey',
                        orientation='horizontal', label='Prior')
        if len(post_filtered) > 0:
            ax_hist.hist(post_filtered, bins=30, alpha=0.5, color='blue',
                        orientation='horizontal', label='Posterior')

        ax_hist.set_xlabel('Count')
        ax_hist.grid(alpha=0.3)
        plt.setp(ax_hist.get_yticklabels(), visible=False)

        if idx == 0:
            ax_hist.legend(fontsize=8)

    plt.tight_layout()

    output_file = os.path.join(OUTPUT_DIR, 'timeseries_heads_prior_posterior.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")


def plot_spring_flux():
    """
    Plot spring flux time series comparing prior and posterior with histograms.

    Data structure: rows = realizations, columns = observations
    """
    print("\n5. Plotting spring flux time series...")

    # Load PST to get headers
    pst_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.pst')
    pst = pyemu.Pst(pst_file)
    headers = pst.observation_data.index.tolist()

    # Find prior and posterior files
    prior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.0.obs.csv')
    par_files = sorted(Path(MASTER_IES_DIR).glob(f'{MODEL_NAME}.*.obs.csv'))
    iter_numbers = [int(f.stem.split('.')[-2]) for f in par_files if f.stem.split('.')[-2].isdigit()]
    last_iter = max(iter_numbers)
    posterior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.{last_iter}.obs.csv')

    # Look for flux observations
    col_startswith = 'ts-flux'
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("  No spring flux observations found, skipping")
        return

    print(f"  Loading {len(desired_cols)} flux observations...")

    # Read data - rows are realizations, columns are observations
    pr_oe_ts = pd.read_csv(prior_file, usecols=desired_cols, index_col=0)
    pt_oe_ts = pd.read_csv(posterior_file, usecols=desired_cols, index_col=0)

    # Extract metadata from column names
    col_metadata = []
    for col in desired_cols:
        time_val = extract_time(col)
        usecol_val = extract_usecol(col)
        col_metadata.append({'col': col, 'time': time_val, 'usecol': usecol_val})

    metadata_df = pd.DataFrame(col_metadata)

    # Get unique locations
    unique_usecols = metadata_df['usecol'].unique()

    if len(unique_usecols) == 0:
        print("  No flux locations found, skipping")
        return

    usecol = unique_usecols[0]  # Use first location
    print(f"  Plotting flux for location: {usecol}")

    # Get columns for this location
    usecol_cols = metadata_df[metadata_df['usecol'] == usecol].sort_values('time')
    time_values = usecol_cols['time'].values
    obs_cols = usecol_cols['col'].values

    # Get data for this location from both prior and posterior
    prior_data = pr_oe_ts[obs_cols].replace(-999, np.nan)  # shape: (n_realizations, n_times)
    post_data = pt_oe_ts[obs_cols].replace(-999, np.nan)   # shape: (n_realizations, n_times)

    # Calculate percentiles for posterior (1st to 99th) across all times and realizations
    all_post_values = post_data.values.flatten()
    all_post_values = all_post_values[~np.isnan(all_post_values)]

    if len(all_post_values) > 0:
        y_min = np.percentile(all_post_values, 1)
        y_max = np.percentile(all_post_values, 99)
    else:
        y_min, y_max = None, None

    # Create figure with gridspec (time series + histogram)
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], figure=fig)

    ax_ts = fig.add_subplot(gs[0, 0])
    ax_hist = fig.add_subplot(gs[0, 1], sharey=ax_ts)

    # Plot prior ensemble (grey) - each row is a realization
    for real_idx in prior_data.index:
        y_vals = prior_data.loc[real_idx].values
        # Only plot if within percentile range
        if y_min is not None and y_max is not None:
            if np.any((y_vals >= y_min) & (y_vals <= y_max) & ~np.isnan(y_vals)):
                ax_ts.plot(time_values, y_vals, alpha=0.2, color='grey', linewidth=0.5)
        else:
            ax_ts.plot(time_values, y_vals, alpha=0.2, color='grey', linewidth=0.5)

    # Plot posterior ensemble (blue) - each row is a realization
    for real_idx in post_data.index:
        y_vals = post_data.loc[real_idx].values
        # Only plot if within percentile range
        if y_min is not None and y_max is not None:
            if np.any((y_vals >= y_min) & (y_vals <= y_max) & ~np.isnan(y_vals)):
                ax_ts.plot(time_values, y_vals, alpha=0.3, color='blue', linewidth=0.8)
        else:
            ax_ts.plot(time_values, y_vals, alpha=0.3, color='blue', linewidth=0.8)

    # Set y-limits based on posterior 1-99 percentiles
    if y_min is not None and y_max is not None:
        y_range = y_max - y_min
        ax_ts.set_ylim(y_min - 0.1 * y_range, y_max + 0.1 * y_range)

    ax_ts.set_title(f'Spring Flux - {usecol}', fontsize=12, fontweight='bold')
    ax_ts.set_ylabel('Flux (m³/d)')
    ax_ts.set_xlabel('Time')
    ax_ts.grid(alpha=0.3)
    ax_ts.legend(['Prior', 'Posterior'], fontsize=8)

    # RIGHT PLOT: Histogram
    # Collect all values for histograms (all times, all realizations)
    prior_all_vals = prior_data.values.flatten()
    prior_all_vals = prior_all_vals[~np.isnan(prior_all_vals)]

    post_all_vals = post_data.values.flatten()
    post_all_vals = post_all_vals[~np.isnan(post_all_vals)]

    # Filter to 1-99 percentile range
    if y_min is not None and y_max is not None:
        prior_filtered = prior_all_vals[(prior_all_vals >= y_min) & (prior_all_vals <= y_max)]
        post_filtered = post_all_vals[(post_all_vals >= y_min) & (post_all_vals <= y_max)]
    else:
        prior_filtered = prior_all_vals
        post_filtered = post_all_vals

    # Plot histograms horizontally
    if len(prior_filtered) > 0:
        ax_hist.hist(prior_filtered, bins=30, alpha=0.5, color='grey',
                    orientation='horizontal', label='Prior')
    if len(post_filtered) > 0:
        ax_hist.hist(post_filtered, bins=30, alpha=0.5, color='blue',
                    orientation='horizontal', label='Posterior')

    ax_hist.set_xlabel('Count')
    ax_hist.grid(alpha=0.3)
    plt.setp(ax_hist.get_yticklabels(), visible=False)
    ax_hist.legend(fontsize=8)

    plt.tight_layout()

    output_file = os.path.join(OUTPUT_DIR, 'spring_flux_prior_posterior.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")


def plot_budgets():
    """
    Plot budget components as boxplots grouped by kper.
    """
    print("\n6. Plotting budget components...")

    # Load PST to get headers
    pst_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.pst')
    pst = pyemu.Pst(pst_file)
    headers = pst.observation_data.index.tolist()

    # Find posterior file
    par_files = sorted(Path(MASTER_IES_DIR).glob(f'{MODEL_NAME}.*.obs.csv'))
    iter_numbers = [int(f.stem.split('.')[-2]) for f in par_files if f.stem.split('.')[-2].isdigit()]
    last_iter = max(iter_numbers)
    posterior_file = os.path.join(MASTER_IES_DIR, f'{MODEL_NAME}.{last_iter}.obs.csv')

    # Look for budget observations
    col_startswith = 'budget'
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("  No budget observations found, skipping")
        return

    print(f"  Loading {len(desired_cols)} budget observations...")

    # Read budget data
    pt_oe_budget = pd.read_csv(posterior_file, usecols=desired_cols)

    # Format the dataframe
    budget_df = pt_oe_budget.transpose()
    budget_df['usecol'] = budget_df.index.map(extract_usecol)
    budget_df['kper'] = budget_df.index.map(extract_kper)

    # Filter out any rows with missing kper
    budget_df = budget_df[budget_df['kper'].notna()]

    unique_usecols = sorted(budget_df['usecol'].unique())
    unique_kpers = sorted(budget_df['kper'].unique())

    print(f"  Found {len(unique_usecols)} budget components across {len(unique_kpers)} stress periods")

    if len(unique_usecols) == 0:
        print("  No budget components found, skipping plot")
        return

    # Create subplots
    n_usecols = len(unique_usecols)
    ncols = min(3, n_usecols)
    nrows = int(np.ceil(n_usecols / ncols))

    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    fig.suptitle('Budget Components by Stress Period (Posterior)', fontsize=14, fontweight='bold')

    if n_usecols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for idx, usecol in enumerate(unique_usecols):
        ax = axes[idx]

        # Get data for this budget component
        usecol_data = budget_df[budget_df['usecol'] == usecol]

        # Prepare data for boxplot (one box per kper)
        box_data = []
        box_labels = []

        for kper in unique_kpers:
            kper_data = usecol_data[usecol_data['kper'] == kper]
            if len(kper_data) > 0:
                # Get value columns (exclude usecol and kper)
                value_cols = [col for col in kper_data.columns if col not in ['usecol', 'kper']]
                values = kper_data[value_cols].values.flatten()
                # Remove NaN values
                values = values[~np.isnan(values)]
                if len(values) > 0:
                    box_data.append(values)
                    box_labels.append(f'kper {int(kper)}')

        if len(box_data) > 0:
            bp = ax.boxplot(box_data, tick_labels=box_labels, patch_artist=True)
            for patch in bp['boxes']:
                patch.set_facecolor('lightblue')
                patch.set_alpha(0.7)

        ax.set_title(usecol, fontsize=10)
        ax.set_ylabel('Value')
        ax.grid(alpha=0.3, axis='y')
        ax.tick_params(axis='x', rotation=45)

    # Hide unused subplots
    for idx in range(n_usecols, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()

    output_file = os.path.join(OUTPUT_DIR, 'budgets_by_kper.png')
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")


def main():
    """Main execution function."""
    print("\n" + "="*70)
    print(f"POST-PROCESSING: {MODEL_NAME}")
    print("="*70 + "\n")

    # Check if master_ies directory exists
    if not os.path.exists(MASTER_IES_DIR):
        print(f"ERROR: Master IES directory not found: {MASTER_IES_DIR}")
        return

    print(f"Master IES directory: {MASTER_IES_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")

    # 1. Plot phi progression
    plot_phi_progression()

    # 2. Plot observation ensembles
    plot_observation_ensembles()

    # 3. Plot parameter ensembles
    plot_parameter_ensembles()

    # 4. Plot time series observations
    plot_timeseries_observations()

    # 5. Plot spring flux
    plot_spring_flux()

    # 6. Plot budgets
    plot_budgets()

    print("\n" + "="*70)
    print("POST-PROCESSING COMPLETE")
    print("="*70 + "\n")


if __name__ == '__main__':
    main()
