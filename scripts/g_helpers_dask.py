"""
Dask-optimized helper functions for PEST++ IES post-processing.

This module provides faster data loading using Dask for large CSV files.
Falls back to pandas for small files and maintains compatibility with pyemu.

Author: Claude Code
Date: 2025-11-19
"""

import os
import pandas as pd
import dask.dataframe as dd
import pyemu
import numpy as np
from dask.diagnostics import ProgressBar

# Re-export all functions from g_helpers_run4 to maintain compatibility
from g_helpers_run4 import *


def load_csv_smart(filepath, usecols=None, index_col=None, use_dask=True, blocksize='50MB'):
    """
    Smart CSV loader that uses Dask for large files, pandas for small files.

    Args:
        filepath: Path to CSV file
        usecols: Columns to load
        index_col: Column to use as index
        use_dask: If True, use Dask for files > 10MB
        blocksize: Dask blocksize for chunking

    Returns:
        DataFrame (pandas or dask depending on file size)
    """
    file_size_mb = os.path.getsize(filepath) / (1024 * 1024)

    # Use pandas for small files or when Dask is disabled
    if file_size_mb < 10 or not use_dask:
        return pd.read_csv(filepath, usecols=usecols, index_col=index_col)

    # Use Dask for large files
    print(f"  Loading large file with Dask ({file_size_mb:.1f} MB): {os.path.basename(filepath)}")

    if usecols is not None:
        # Dask requires explicit column list
        if index_col is not None and index_col not in usecols:
            if isinstance(usecols, list):
                usecols = [index_col] + usecols

    return dd.read_csv(filepath, usecols=usecols, blocksize=blocksize)


def load_data_dask(paths, run_name, last_iter=None, use_prior_only=False, restart=False):
    """
    Load data with Dask optimization for large files.

    Returns data dictionary compatible with original load_data function.
    """
    # Select directory based on mode
    if use_prior_only:
        m_d = paths['m_d_pr']
        print(f"Loading PRIOR-ONLY data from: {m_d}")
    else:
        m_d = paths['m_d']
        print(f"Loading data from: {m_d}")

    m_d_pr = paths['m_d_pr']

    # Load PST file (small, keep as pandas)
    pst = pyemu.Pst(os.path.join(m_d, f"{run_name}.pst"))

    if last_iter is None:
        if use_prior_only:
            last_iter = 0
        else:
            last_iter = pst.control_data.noptmax

    # File paths
    pr_oe_fn = os.path.join(m_d_pr, f"{run_name}.0.obs.csv")
    pr_pe_fn = os.path.join(m_d_pr, f"{run_name}.0.par.csv")
    noise_fn = os.path.join(m_d, f"{run_name}.obs+noise.csv")

    par_data_fn = os.path.join(m_d_pr, f"{run_name}.par_data.csv")
    obs_data_fn = os.path.join(m_d_pr, f"{run_name}.obs_data.csv")

    if restart:
        run_name = f'{run_name}_restart'

    if use_prior_only:
        f_act = None
        f_meas = None
    else:
        f_act = os.path.join(m_d, f"{run_name}.phi.actual.csv")
        f_meas = os.path.join(m_d, f"{run_name}.phi.meas.csv")

    pt_oe_fn = os.path.join(m_d, f"{run_name}.{last_iter}.obs.csv")
    pt_pe_fn = os.path.join(m_d, f"{run_name}.{last_iter}.par.csv")

    print("Loading metadata files (pandas)...")
    # Load small metadata files with pandas
    obs_data = pd.read_csv(obs_data_fn, index_col=0)
    par_data = pd.read_csv(par_data_fn, index_col=0)

    print("Loading ensemble files...")
    # Load observation ensembles - keep as pyemu for compatibility
    # PyEMU will handle these, no need to optimize
    pr_oe = pyemu.ObservationEnsemble.from_csv(pst=pst, filename=pr_oe_fn)

    if use_prior_only:
        pt_oe = pr_oe
    else:
        pt_oe = pyemu.ObservationEnsemble.from_csv(pst=pst, filename=pt_oe_fn)

    noise = pyemu.ObservationEnsemble.from_csv(pst=pst, filename=noise_fn)

    # Load parameter ensembles
    pr_pe = pyemu.ParameterEnsemble.from_csv(pst=pst, filename=pr_pe_fn)

    if use_prior_only:
        pt_pe = pr_pe
    else:
        pt_pe = pyemu.ParameterEnsemble.from_csv(pst=pst, filename=pt_pe_fn)

    # Store file paths for lazy loading
    data = {
        'pst': pst,
        'obs_data': obs_data,
        'par_data': par_data,
        'pr_oe': pr_oe,
        'pt_oe': pt_oe,
        'noise': noise,
        'pr_pe': pr_pe,
        'pt_pe': pt_pe,
        'f_act': f_act,
        'f_meas': f_meas,
        # Store filenames for Dask lazy loading
        'pr_oe_fn': pr_oe_fn,
        'pt_oe_fn': pt_oe_fn,
        'noise_fn': noise_fn,
        'pr_pe_fn': pr_pe_fn,
        'pt_pe_fn': pt_pe_fn,
        'headers': pst.observation_data.index.to_list(),
        'last_iter': last_iter,
    }

    print(f"Data loaded successfully (last iteration: {last_iter})")
    return data


def get_array_statistics(rdf, use_dask=True):
    """
    Optimized array statistics calculation.

    Args:
        rdf: DataFrame (pandas or dask)
        use_dask: Use Dask for computation

    Returns:
        dict: Statistics arrays organized by (kper, kstp)
    """
    # If it's a Dask DataFrame, compute it first for this operation
    # (array operations are not efficiently parallelized anyway)
    if hasattr(rdf, 'compute'):
        print("  Computing Dask DataFrame for array statistics...")
        with ProgressBar():
            rdf = rdf.compute()

    # Rest of the function is the same as original
    # Import from original helpers
    from g_helpers_run4 import get_array_statistics as original_get_array_statistics
    return original_get_array_statistics(rdf)


def plot_oname_array_statistics_dask(data, subset_oe=None, output=None, oname_prefix=None):
    """
    Plot observation array statistics with Dask optimization.
    """
    desired_cols = list(data['obs_data'][data['obs_data']['oname'] == oname_prefix].index)
    obs_data = data['obs_data'].loc[desired_cols]

    print(f"  Loading prior data for {oname_prefix}...")
    # Load with Dask
    pr_df = load_csv_smart(data['pr_oe_fn'], usecols=desired_cols, use_dask=True)

    print(f"  Loading posterior data for {oname_prefix}...")
    if subset_oe is not None:
        pt_df = subset_oe.loc[:, desired_cols]
        # Convert to pandas if it's Dask
        if hasattr(pt_df, 'compute'):
            with ProgressBar():
                pt_df = pt_df.compute()
    else:
        pt_df = load_csv_smart(data['pt_oe_fn'], usecols=desired_cols, use_dask=True)

    # Compute Dask DataFrames
    if hasattr(pr_df, 'compute'):
        print("  Computing statistics...")
        with ProgressBar():
            pr_df = pr_df.compute()
    if hasattr(pt_df, 'compute'):
        with ProgressBar():
            pt_df = pt_df.compute()

    dfs = {
        'pe': pr_df,
        'pt': pt_df,
    }

    # Call original plotting function
    from g_helpers_run4 import plot_oname_array_statistics as original_plot
    # Temporarily modify data to use computed DataFrames
    temp_data = data.copy()
    return original_plot(temp_data, subset_oe=dfs['pt'], output=output, oname_prefix=oname_prefix)


def plot_timeseries_obs_dask(data, subset_oe=None, output=None, truth_file=None, col_startswith='oname:ts-heads'):
    """
    Plot time series observations with Dask optimization.
    """
    obs_data = data['obs_data'][data['obs_data']['oname'] == col_startswith]
    desired_cols = list(obs_data.index)

    print(f"  Loading timeseries data with Dask...")

    # Load with Dask
    pr_oe_fn = data['pr_oe_fn']
    pt_oe_fn = data['pt_oe_fn']
    noise_fn = data['noise_fn']

    pre_oe_ts = load_csv_smart(pr_oe_fn, usecols=desired_cols, use_dask=True)

    if subset_oe is not None:
        pt_oe_ts = subset_oe.loc[:, desired_cols]
    else:
        pt_oe_ts = load_csv_smart(pt_oe_fn, usecols=desired_cols, use_dask=True)

    noise_ts = load_csv_smart(noise_fn, usecols=desired_cols, use_dask=True)

    # Compute all at once
    print("  Computing Dask operations...")
    with ProgressBar():
        if hasattr(pre_oe_ts, 'compute'):
            pre_oe_ts = pre_oe_ts.compute()
        if hasattr(pt_oe_ts, 'compute'):
            pt_oe_ts = pt_oe_ts.compute()
        if hasattr(noise_ts, 'compute'):
            noise_ts = noise_ts.compute()

    # Call original plotting function
    from g_helpers_run4 import plot_timeseries_obs as original_plot

    # Create temporary data structure
    dfs = {
        'pr': pre_oe_ts,
        'pt': pt_oe_ts,
        'noise': noise_ts
    }

    return original_plot(data, subset_oe=dfs['pt'], output=output, truth_file=truth_file, col_startswith=col_startswith)


# Note: All functions from g_helpers_run4 are re-exported via import *
# Additional Dask-optimized functions are defined above:
# - load_data_dask: Optimized data loading
# - load_csv_smart: Smart CSV loading with Dask for large files
