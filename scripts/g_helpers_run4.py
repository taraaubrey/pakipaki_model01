import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from scipy import stats
from scipy.interpolate import griddata
import pyemu
import flopy

from setup import MODEL_NAME, NREALS_PRIOR

plt.rcParams['font.family'] = 'Arial'  # or 'Times New Roman', 'Helvetica', etc.
plt.rcParams['font.size'] = 8

# Extraction helper functions
extract_kper = lambda x: int([s for s in x.split('_') if s.startswith('kper:')][0].split(':')[1])
extract_kstp = lambda x: int([s for s in x.split('_') if s.startswith('kstp:')][0].split(':')[1])

extract_hkper = lambda x: int([s for s in x.split('-') if s.startswith('kper')][0][4:])
extract_hkstp = lambda x: int([s for s in x.split('-') if s.startswith('kstp')][0].split('_')[0][4:])

extract_i = lambda x: int([s for s in x.split('_') if s.startswith('i:')][0].split(':')[1])
extract_j = lambda x: int([s for s in x.split('_') if s.startswith('j:')][0].split(':')[1])
extract_idx1 = lambda x: int([s for s in x.split('_') if s.startswith('idx1:')][0].split(':')[1])
extract_idx2 = lambda x: int([s for s in x.split('_') if s.startswith('idx2:')][0].split(':')[1])
extract_time = lambda x: int([s for s in x.split('_') if s.startswith('time:')][0].split(':')[1])
extract_usecol = lambda x: [s for s in x.split('_') if s.startswith('usecol:')][0].split(':')[1]

def setup_paths(run_name, m_d_str='master_ies', m_d_prior=None):
    """Set up directory paths for a given run."""
    # Get absolute path from scripts directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)

    paths = {
        'o_d': os.path.join(project_root, f'models/{run_name}/{run_name}'),
        'temp_d': os.path.join(project_root, f'models/{run_name}/{run_name}'),
        'm_d': os.path.join(project_root, f'models/{run_name}/pest/{m_d_str}'),
        'output': os.path.join(project_root, f'models/{run_name}/figures'),
    }

    if m_d_prior:
        paths['m_d_pr'] = os.path.join(project_root, f'models/{run_name}/pest/{m_d_prior}')
    else:
        paths['m_d_pr'] = os.path.join(project_root, f'models/{run_name}/pest/{m_d_str}')
    # Create output directory if it doesn't exist
    os.makedirs(paths['output'], exist_ok=True)

    return paths

def load_data(paths, run_name, last_iter=None, use_prior_only=False, restart=False):
    """Load all necessary data files.

    Args:
        paths: Dictionary of paths
        run_name: Name of the model run
        last_iter: Last iteration number (for history matching runs)
        use_prior_only: If True, load only prior data from master_ies_prior directory
    """
    # Select directory based on mode
    if use_prior_only:
        m_d = paths['m_d_pr']
        print(f"Loading PRIOR-ONLY data from: {m_d}")
    else:
        m_d = paths['m_d']
        print(f"Loading data from: {m_d}")
    
    m_d_pr = paths['m_d_pr']
    
    # Load PST file
    pst = pyemu.Pst(os.path.join(m_d, f"{run_name}.pst"))

    if last_iter is None:
        if use_prior_only:
            last_iter = 0  # Prior only has iteration 0
        else:
            last_iter = pst.control_data.noptmax

    # prior first
    pr_oe_fn = os.path.join(m_d_pr, f"{run_name}.0.obs.csv")  # observations
    pr_pe_fn = os.path.join(m_d_pr, f"{run_name}.0.par.csv")  # parameters
    noise_fn = os.path.join(m_d, f"{run_name}.obs+noise.csv") # obs + noise

    par_data_fn = os.path.join(m_d_pr, f"{run_name}.par_data.csv")
    obs_data_fn = os.path.join(m_d_pr, f"{run_name}.obs_data.csv") 

    if restart:
        run_name = f'{MODEL_NAME}_restart'

    # File paths (only available for history matching runs)
    if use_prior_only:
        f_act = None
        f_meas = None
    else:
        f_act = os.path.join(m_d, f"{run_name}.phi.actual.csv")
        f_meas = os.path.join(m_d, f"{run_name}.phi.meas.csv")

    # Posterior
    pt_oe_fn = os.path.join(m_d, f"{run_name}.{last_iter}.obs.csv") # observations
    pt_pe_fn = os.path.join(m_d, f"{run_name}.{last_iter}.par.csv") # parameters

    # # Load observation data
    obs_data = pd.read_csv(obs_data_fn, index_col=0)
    par_data = pd.read_csv(par_data_fn, index_col=0)

    # # Load observation ensembles
    pr_oe = pyemu.ObservationEnsemble.from_csv(pst=pst, filename=pr_oe_fn)

    # For prior-only mode, pt_oe is same as pr_oe
    # if use_prior_only:
    #     pt_oe = pr_oe
    # else:
    #     pt_oe = pyemu.ObservationEnsemble.from_csv(pst=pst, filename=pt_oe_fn)

    noise = pyemu.ObservationEnsemble.from_csv(pst=pst, filename=noise_fn)

    # Load parameter ensembles
    pr_pe = pyemu.ParameterEnsemble.from_csv(pst=pst, filename=pr_pe_fn)

    # For prior-only mode, pt_pe is same as pr_pe
    # if use_prior_only:
    #     pt_pe = pr_pe
    # else:
    #     pt_pe = pyemu.ParameterEnsemble.from_csv(pst=pst, filename=pt_pe_fn)

    # Get headers
    headers = obs_data.index.to_list()
    param_headers = par_data.index.tolist()

    data = {
        'pst': pst,
        'f_act': f_act,
        'f_meas': f_meas,
        'pr_oe_fn': pr_oe_fn,
        'pt_oe_fn': pt_oe_fn,
        'pr_pe_fn': pr_pe_fn,
        'pt_pe_fn': pt_pe_fn,
        'noise_fn': noise_fn,
        'obs_data': obs_data,
        'par_data': par_data,
        'pr_oe': pr_oe,
        # 'pt_oe': pt_oe,
        'noise': noise,
        'pr_pe': pr_pe,
        # 'pt_pe': pt_pe,
        'headers': headers,
        'param_headers': param_headers,
        'use_prior_only': use_prior_only,
    }

    return data

def subset_realizations_by_phi(data, min_threshold=100):
    phi = pd.read_csv(data['f_act'], index_col=0)

    # Convert to numpy array, flatten, get argsort, then convert back to DataFrame indexes
    flat_values = phi.iloc[:,5:-1].values.flatten()
    flat_indexes = np.argsort(flat_values)[:min_threshold]

    # Convert flat indexes back to DataFrame row/column positions
    rows, cols = np.unravel_index(flat_indexes, phi.iloc[:,5:-1].shape)
    int_rows = [int(row) for row in rows] # iterations
    int_cols = [int(col) for col in cols] # realizations
    indexes = list(zip(int_rows, int_cols))
    vals = [int(phi.iloc[:,5:-1].iloc[r, c]) for r, c in indexes]
    return indexes, vals

def get_subset_ensembles(paths, run_name, subset_indexes):
    m_d = paths['m_d']
    unique_iterations = set(i[0] for i in subset_indexes)

    oe_list=[]
    pe_list=[]
    for iter in unique_iterations:
        oe_fn = os.path.join(m_d, f"{run_name}.{iter}.obs.csv")
        pe_fn = os.path.join(m_d, f"{run_name}.{iter}.par.csv")
        # get the realization indexes
        reals = [i[1] for i in subset_indexes if i[0] == iter]
        obs = pd.read_csv(
            oe_fn,
            index_col=0)
        pars = pd.read_csv(
            pe_fn,
            index_col=0)
        
        #convert index to string
        reals = [str(r) for r in reals]
        obs.index = obs.index.map(str)
        pars.index = pars.index.map(str)

        oe_list.append(obs.loc[reals,:])
        pe_list.append(pars.loc[reals,:])
    
    # merge all dataframes
    subset_oe = pd.concat(oe_list)
    subset_pe = pd.concat(pe_list)

    # save
    subset_oe.to_csv(os.path.join(paths['output'], 'subset_oe.csv'))
    subset_pe.to_csv(os.path.join(paths['output'], 'subset_pe.csv'))

    return subset_oe, subset_pe




def plot_phi_comparison(data, output):
    """Plot actual vs measured phi over iterations.

    Note: Only available for history matching runs (not prior-only mode).
    """
    if data.get('use_prior_only', False):
        print("Skipping phi comparison plot (not available for prior-only mode)")
        return

    f_act = data['f_act']
    f_meas = data['f_meas']

    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10, 3.5))

    # Left - Actual phi
    ax = axes[0]
    phi = pd.read_csv(f_act, index_col=0)
    phi.index = phi.total_runs
    phi.iloc[:, 5:].apply(np.log10).plot(legend=False, lw=0.5, color='k', ax=ax)
    ax.set_title(r'Actual $\Phi$')
    ax.set_ylabel(r'log $\Phi$')
    ax.set_xlabel('Total Model Runs')

    # Right - Measured phi
    ax = axes[-1]
    phi = pd.read_csv(f_meas, index_col=0)
    phi.index = phi.total_runs
    phi.iloc[:, 6:].apply(np.log10).plot(legend=False, lw=0.2, color='r', ax=ax)
    ax.set_title(r'Measured+Noise $\Phi$')
    ax.set_xlabel('Total Model Runs')

    fig.tight_layout()
    fig_fn = os.path.join(output, 'output.phi_comparison.png')
    plt.savefig(fig_fn, dpi=300)
    plt.close()
    print(f"Saved: {fig_fn}")

def plot_phi_boxplot(data, output):
    """Plot boxplot of phi distribution by iteration.

    Note: Only available for history matching runs (not prior-only mode).
    """
    if data.get('use_prior_only', False):
        print("Skipping phi boxplot (not available for prior-only mode)")
        return

    f_act = data['f_act']
    phi_act = pd.read_csv(f_act, index_col=0)

    # Prepare data for boxplot
    box_data = []
    value_columns = phi_act.columns[5:].to_list()

    for total_runs in phi_act.index:
        row_values = []
        for col in value_columns:
            val = phi_act.loc[total_runs, col]
            if pd.notna(val):
                row_values.append(val)
        box_data.append(row_values)

    # Create boxplot
    fig, ax = plt.subplots(figsize=(10, 6))
    bp = ax.boxplot(box_data, patch_artist=True)

    # Color boxes
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
        patch.set_alpha(0.7)

    ax.set_xlabel('IES Iterations')
    ax.set_ylabel('Phi Values')
    ax.set_title('Distribution of Phi Values by Total Runs')
    ax.set_yscale('log')
    ax.grid(which='both', alpha=0.3)

    # Set custom x-axis labels (iterations)
    ax.set_xticklabels([str(i) for i in range(len(box_data))])

    plt.tight_layout()
    fig_fn = os.path.join(output, 'output.phi_by_total_runs_boxplot.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

def format_tsdf(df):
    """Format time series dataframe."""
    tdf = df.transpose()
    tdf['time'] = tdf.index.map(extract_time)
    tdf['kper'] = tdf.index.map(extract_kper)
    tdf['kstp'] = tdf.index.map(extract_kstp)
    tdf['usecol'] = tdf.index.map(extract_usecol)
    return tdf

def plot_timeseries_obs(data, subset_oe, output, truth_file=None, col_startswith='oname:ts-heads'):
    """Plot time series observations - prior vs posterior vs truth.

    Creates two separate plots:
    1. obs.timeseries_heads.png - times 1-53 with truth values (all ensembles)
    2. prediction.timeseries_heads.png - times 54+ without truth values (all ensembles)
    """
    headers = data['headers']
    pr_oe_fn = data['pr_oe_fn']
    noise_fn = data['noise_fn']

    # Load truth data
    if truth_file is None:
        truth_df = None
    else:
        if os.path.exists(truth_file):
            truth_df = pd.read_csv(truth_file, index_col=0)
        else:
            print(f"Warning: Truth file not found: {truth_file}")
            truth_df = None

    # Get time series head observations
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("No time series head observations found, skipping plot")
        return

    pre_oe_ts = pd.read_csv(pr_oe_fn, usecols=desired_cols)
    pt_oe_ts = subset_oe.loc[:, list(desired_cols)]
    noise_ts = pd.read_csv(noise_fn, usecols=desired_cols)

    ts_df = format_tsdf(pre_oe_ts)
    ts_pt_df = format_tsdf(pt_oe_ts)
    ts_noise_df = format_tsdf(noise_ts)

    unique_usecols = ts_df['usecol'].unique()
    n = len(unique_usecols)

    # ========== PLOT 1: Times 1-53 with truth values ==========
    fig1 = plt.figure(figsize=(14, n * 3))
    gs1 = gridspec.GridSpec(n, 2, width_ratios=[4, 1])

    axes1 = []
    for i in range(n):
        if i == 0:
            ax0 = fig1.add_subplot(gs1[i, 0])
            ax1 = fig1.add_subplot(gs1[i, 1], sharey=ax0)
        else:
            ax0 = fig1.add_subplot(gs1[i, 0], sharex=axes1[0][0])
            ax1 = fig1.add_subplot(gs1[i, 1], sharey=ax0)
        axes1.append([ax0, ax1])
        

    for idx, usecol in enumerate(unique_usecols):
        
        # Get truth if available
        if truth_df is not None and usecol in truth_df.columns:
            truth_col = truth_df[usecol]
        else:
            truth_col = None

        # Prep data - filter for times 1-53
        mask = (ts_df['usecol'] == usecol) & (ts_df['kper'] == 2)
        df_pr_kper2 = ts_df[mask].sort_values('kstp').replace(-999, np.nan)
        df_pt_kper2 = ts_pt_df[mask].sort_values('kstp').replace(-999, np.nan)
        df_obsnoise_kper2 = ts_noise_df[mask].sort_values('kstp').replace(-999, np.nan)

        value_cols = [col for col in df_pt_kper2.columns if col not in ['time', 'kper', 'kstp', 'usecol']]

        # Calculate y-limits
        pt_values = df_pt_kper2[value_cols].values.flatten()
        if len(pt_values) > 0:
            y_min = np.min(pt_values)
            y_max = np.max(pt_values)
        else:
            y_min, y_max = None, None

        # LEFT PLOT: Time series (NO outlier removal - plot all ensembles)
        ax_ts = axes1[idx][0]

        # Plot all ensemble members
        for col in value_cols:
            y_pre = df_pr_kper2[col].copy() if col in df_pr_kper2.columns else None
            if y_pre is not None:
                ax_ts.plot(df_pr_kper2['kstp'], y_pre, alpha=0.2, color='grey', linewidth=0.5)
            y_pt = df_pt_kper2[col].copy() if col in df_pt_kper2.columns else None
            if y_pt is not None:
                ax_ts.plot(df_pt_kper2['kstp'], y_pt, alpha=0.3, color='blue', linewidth=0.8)

        # Plot truth
        if truth_col is not None:
            truth_col_filtered = truth_col.iloc[:53]  # Only times 1-53
            ax_ts.plot(truth_col_filtered, alpha=1, color='red', linewidth=2, label='Truth')

        ax_ts.set_title(f'{usecol}', loc='left')
        ax_ts.set_ylabel('Head (m)')
        ax_ts.grid(alpha=0.3)

        # Set y-limits based on prior 5-95th percentiles
        if y_min is not None and y_max is not None:
            ax_ts.set_ylim(y_min, y_max)

        if idx == 0 and truth_col is not None:
            ax_ts.legend(['Prior', 'Posterior', 'Truth'], fontsize=8)

        # Only show x-label on bottom plot
        if idx == n - 1:
            ax_ts.set_xlabel('Days')
        else:
            ax_ts.set_xlabel('')


        # RIGHT PLOT: Histogram for kper 1-2 (present) - ONLY ONE HISTOGRAM
        df_pr_kper2_vals = df_pr_kper2[value_cols].values.flatten()
        df_pt_kper2_vals = df_pt_kper2[value_cols].values.flatten()
        df_obsnoise_kper2_vals = df_obsnoise_kper2[value_cols].values.flatten()
        # df_pr_kper4_vals = df_pr_kper4[value_cols].values.flatten()
        # df_pt_kper4_vals = df_pt_kper4[value_cols].values.flatten()
        
        ax_hist_kper12 = axes1[idx][1]

        if len(df_pr_kper2_vals) > 0:
            ax_hist_kper12.hist(df_pr_kper2_vals, bins=20, alpha=0.5, color='lightgrey',
                        orientation='horizontal', density=True, edgecolor='grey', linewidth=0.8,
                        label='Prior')

        if len(df_pt_kper2_vals) > 0:
            ax_hist_kper12.hist(df_pt_kper2_vals, bins=20, alpha=0.6, color='lightblue',
                        orientation='horizontal', density=True, edgecolor='blue', linewidth=0.8,
                        label='Posterior')

        # Add obs+noise histogram if available
        if len(df_obsnoise_kper2_vals) > 0:
            ax_hist_kper12.hist(
                df_obsnoise_kper2_vals, bins=20, alpha=0.4, color='orange',
                orientation='horizontal', density=True, edgecolor='darkorange', linewidth=0.8,
                label='Obs+noise')

        # Add mean line for "present" (kper 1-2)
        if len(df_obsnoise_kper2_vals) > 0:
            mean_pt = np.mean(df_pt_kper2_vals)
            mean_truth = np.mean(df_obsnoise_kper2_vals)
            ax_hist_kper12.axhline(mean_pt, color='green', linestyle='-', linewidth=2, label='Mean')
            ax_hist_kper12.axhline(mean_truth, color='red', linestyle='--', linewidth=2, label='Obs+Noise Mean')

        ax_hist_kper12.set_xlabel('Density', fontsize=8)
        ax_hist_kper12.set_title('kper 2', fontsize=9)
        ax_hist_kper12.tick_params(axis='both', labelsize=8)
        ax_hist_kper12.grid(alpha=0.3)
        plt.setp(ax_hist_kper12.get_yticklabels(), visible=False)

        # Add legend only to first subplot
        if idx == 0:
            ax_hist_kper12.legend(fontsize=6, loc='upper right')

    plt.tight_layout()
    fig_fn = os.path.join(output, 'obs.timeseries_heads.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

    # ========== PLOT 2: Times 54+ without truth values ==========
    fig2 = plt.figure(figsize=(18, n * 3))
    gs2 = gridspec.GridSpec(n, 3, width_ratios=[4, 1, 1])

    axes2 = []
    for i in range(n):
        if i == 0:
            ax0 = fig2.add_subplot(gs2[i, 0])
            ax1 = fig2.add_subplot(gs2[i, 1], sharey=ax0)
            ax2 = fig2.add_subplot(gs2[i, 2], sharey=ax0)
        else:
            ax0 = fig2.add_subplot(gs2[i, 0], sharex=axes2[0][0])
            ax1 = fig2.add_subplot(gs2[i, 1], sharey=ax0)
            ax2 = fig2.add_subplot(gs2[i, 2], sharey=ax0)
        axes2.append([ax0, ax1, ax2])

    for idx, usecol in enumerate(unique_usecols):
        mask = (ts_df['usecol'] == usecol) & (ts_df['kper'] == 2)
        df_pr_kper2 = ts_df[mask].sort_values('kstp').replace(-999, np.nan)
        df_pt_kper2 = ts_pt_df[mask].sort_values('kstp').replace(-999, np.nan)
        df_obsnoise_kper2 = ts_noise_df[mask].sort_values('kstp').replace(-999, np.nan)
        # KPER4
        mask = (ts_df['usecol'] == usecol) & (ts_df['kper'] == 4)
        df_pr_kper4 = ts_df[mask].sort_values('kstp').replace(-999, np.nan)
        df_pt_kper4 = ts_pt_df[mask].sort_values('kstp').replace(-999, np.nan)

        value_cols = [col for col in df_pr_kper4.columns if col not in ['kstp', 'usecol']]
        value_cols_pt = [col for col in df_pt_kper4.columns if col not in ['kstp', 'usecol']]

        # LEFT PLOT: Time series (NO outlier removal - plot all ensembles)
        ax_ts = axes2[idx][0]

        # Plot all ensemble members
        for col in value_cols:
            y = df_pr_kper4[col].copy() if col in df_pr_kper4.columns else None
            if y is not None:
                ax_ts.plot(df_pr_kper4['kstp'], y, alpha=0.2, color='grey', linewidth=0.5)

        for col in value_cols_pt:
            y_pt = df_pt_kper4[col].copy() if col in df_pt_kper4.columns else None
            if y_pt is not None:
                ax_ts.plot(df_pt_kper4['kstp'], y_pt, alpha=0.3, color='blue', linewidth=0.8)

        ymin = np.percentile(df_pt_kper4[value_cols_pt].values.flatten(), 1)
        ymax = np.percentile(df_pt_kper4[value_cols_pt].values.flatten(), 99)

        ax_ts.set_title(f'{usecol}', loc='left')
        ax_ts.set_ylabel('Head (m)')
        ax_ts.grid(alpha=0.3)
        ax_ts.set_ylim(ymin, ymax)

        if idx == 0:
            ax_ts.legend(['Prior', 'Posterior'], fontsize=8)

        # Only show x-label on bottom plot
        if idx == n - 1:
            ax_ts.set_xlabel('Time Step')
        else:
            ax_ts.set_xlabel('')
        

        if idx == 0 and truth_col is not None:
            ax_ts.legend(['Prior', 'Posterior', 'Truth'], fontsize=8)

        df_pr_kper2_vals = df_pr_kper2[value_cols].values.flatten()
        df_pt_kper2_vals = df_pt_kper2[value_cols_pt].values.flatten()
        df_obsnoise_kper2_vals = df_obsnoise_kper2[value_cols].values.flatten()
        df_pr_kper4_vals = df_pr_kper4[value_cols].values.flatten()
        df_pt_kper4_vals = df_pt_kper4[value_cols_pt].values.flatten()

        # Determine shared bins for consistent plotting across both histograms
        all_vals = list(df_pr_kper2_vals) + list(df_pt_kper2_vals) + list(df_pr_kper4_vals) + list(df_pt_kper4_vals)
        if len(all_vals) > 0:
            bins = np.linspace(np.percentile(all_vals, 2.5), np.percentile(all_vals, 97.5), 15)
        else:
            bins = 15

        # MIDDLE PLOT: Histogram for first half
        ax_hist_kper12 = axes2[idx][1]

        if len(df_pr_kper2_vals) > 0:
            ax_hist_kper12.hist(df_pr_kper2_vals, bins=bins, alpha=0.5, color='lightgrey',
                        orientation='horizontal', density=True, edgecolor='grey', linewidth=0.8,
                        label='Prior')

        if len(df_pt_kper2_vals) > 0:
            ax_hist_kper12.hist(df_pt_kper2_vals, bins=bins, alpha=0.6, color='lightblue',
                        orientation='horizontal', density=True, edgecolor='blue', linewidth=0.8,
                        label='Posterior')

        # Add mean line
        if len(df_pt_kper2_vals) > 0:
            mean_present = np.mean(df_pt_kper2_vals)
            ax_hist_kper12.axhline(mean_present, color='green', linestyle='-', linewidth=2, label='Mean')

        ax_hist_kper12.set_xlabel('Density', fontsize=8)
        ax_hist_kper12.set_title('First Half', fontsize=9)
        ax_hist_kper12.set_ylim(ymin, ymax)
        ax_hist_kper12.tick_params(axis='both', labelsize=8)
        ax_hist_kper12.grid(alpha=0.3)
        plt.setp(ax_hist_kper12.get_yticklabels(), visible=False)

        # Add legend only to first subplot
        if idx == 0:
            ax_hist_kper12.legend(fontsize=6, loc='upper right')

        # RIGHT PLOT: Histogram for second half
        ax_hist_kper34 = axes2[idx][2]

        if len(df_pr_kper4_vals) > 0:
            ax_hist_kper34.hist(df_pr_kper4_vals, bins=bins, alpha=0.5, color='darkgrey',
                        orientation='horizontal', density=True, edgecolor='black', linewidth=0.8,
                        label='Prior')

        if len(df_pt_kper4_vals) > 0:
            ax_hist_kper34.hist(df_pt_kper4_vals, bins=bins, alpha=0.6, color='darkblue',
                        orientation='horizontal', density=True, edgecolor='navy', linewidth=0.8,
                        label='Posterior')

        # Add mean line
        if len(df_pt_kper4_vals) > 0:
            mean_past = np.mean(df_pt_kper4_vals)
            ax_hist_kper34.axhline(mean_past, color='purple', linestyle='-', linewidth=2, label='Mean')

        ax_hist_kper34.set_xlabel('Density', fontsize=8)
        ax_hist_kper34.set_title('Second Half', fontsize=9)
        ax_hist_kper34.tick_params(axis='both', labelsize=8)
        ax_hist_kper34.grid(alpha=0.3)
        ax_hist_kper34.set_ylim(ymin, ymax)
        plt.setp(ax_hist_kper34.get_yticklabels(), visible=False)

        # Add legend only to first subplot
        if idx == 0:
            ax_hist_kper34.legend(fontsize=7, loc='upper right')

    plt.tight_layout()
    fig_fn = os.path.join(output, 'prediction.timeseries_heads.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")


def plot_budget_histogram(data, output, col_startswith='budget'):
    headers = data['headers']

    # Get budget observations
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("No budget observations found, skipping plot")

    # Read budget data
    budget = {
        'pr': pd.read_csv(data['pr_oe_fn'], usecols=desired_cols),
        'pt':  pd.read_csv(data['pt_oe_fn'], usecols=desired_cols)
    }

    # Format the dataframe to extract usecol and kper
    for name, df in budget.items():
        df = df.transpose()
        df['usecol'] = df.index.map(extract_usecol)
        df['kper'] = df.index.map(extract_kper)
        budget[name] = df
    
    pe_value_cols = [col for col in budget['pr'].columns if col not in ['usecol', 'kper']]
    pt_value_cols = [col for col in budget['pt'].columns if col not in ['usecol', 'kper']]

    for col in budget['pt']['usecol'].unique():
        # create histogram with 4 columns for kper 1-4
        fig, axes = plt.subplots(1, 4, figsize=(7.1, 3))

        for i, ax in enumerate(axes.flatten()):
            mask = (budget['pr']['usecol'] == col) & (budget['pr']['kper'] == i + 1)

            pr_values = budget['pr'][mask][pe_value_cols].values.flatten()
            pt_values = budget['pt'][mask][pt_value_cols].values.flatten()

            if len(pt_values) > 0:
                min = np.percentile(pt_values, 1)
                max = np.percentile(pt_values, 99)
                ax.set_xlim(min, max)

            ax.hist(pr_values, bins=30, alpha=0.5, label='Prior', color='grey')
            ax.hist(pt_values, bins=30, alpha=0.5, label='Posterior', color='blue')
            ax.set_title(f'kper {i + 1}')
        
        axes[0].set_xlabel('Budget Value')
        axes[0].set_ylabel('Frequency')
        axes[0].legend()
        fig.tight_layout()
        fig_fn = os.path.join(output, f'output.budget_{col}.png')
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
        plt.close()
    
    for col in ['awanui', 'confined', 'inflow', 'spring', 'stoss', 'sw']:
        # create histogram with 4 columns for kper 1-4
        fig, axes = plt.subplots(1, 2, figsize=(7.1, 3))

        kpers = [(1, 3), (2, 4)]

        for i, kper_i in enumerate(kpers):
            ax = axes.flatten()[i]
            mask_kper_a = (budget['pt']['usecol'] == col) & (budget['pt']['kper'] == kper_i[0])
            mask_kper_b = (budget['pt']['usecol'] == col) & (budget['pt']['kper'] == kper_i[1])

            kper_a_values = budget['pt'][mask_kper_a][pt_value_cols].values.flatten()
            kper_b_values = budget['pt'][mask_kper_b][pt_value_cols].values.flatten()
            
            # if len(pt_values) > 0:
            #     min = np.percentile(pt_values, 1)
            #     max = np.percentile(pt_values, 99)
            #     ax.set_xlim(min, max)

            ax.hist(kper_a_values, bins=30, alpha=0.5, label='Present', color='green')
            ax.hist(kper_b_values, bins=30, alpha=0.5, label='Past', color='purple')
            ax.set_title(f'kper {i + 1}')
        
        axes[0].set_ylabel('Frequency')
        axes[0].legend()
        fig.tight_layout()
        fig_fn = os.path.join(output, f'prediction.budget_{col}.png')
        plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
        plt.close()


def plot_budgets(data, output, col_startswith='budget'):
    """Plot budget components with subplots for each usecol, showing boxplots by kper.

    Creates a figure with 1 row of subplots (one per usecol), where each subplot
    contains 4 boxplots (one for each kper period).

    Args:
        data: Dictionary containing loaded data
        output: Output directory for saving figure
        col_startswith: Prefix for budget columns (default: 'budget')
    """
    headers = data['headers']
    pt_oe_fn = data['pt_oe_fn']

    # Get budget observations
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("No budget observations found, skipping plot")
        return

    # Read budget data
    pt_oe_budget = pd.read_csv(pt_oe_fn, usecols=desired_cols)

    # Format the dataframe to extract usecol and kper
    budget_df = pt_oe_budget.transpose()
    budget_df['usecol'] = budget_df.index.map(extract_usecol)
    budget_df['kper'] = budget_df.index.map(
        lambda x: int([s for s in x.split('_') if s.startswith('kper:')][0].split(':')[1])if any(s.startswith('kper:') for s in x.split('_')) else 0)

    # Create 'surface water' as sum of awanui, poukawa, and spring
    sw_components = ['awanui', 'poukawa', 'spring']
    sw_data = []

    for kper in sorted(budget_df['kper'].unique()):
        # Get data for each component for this kper
        component_data = []
        for comp in sw_components:
            comp_rows = budget_df[(budget_df['usecol'] == comp) & (budget_df['kper'] == kper)]
            if len(comp_rows) > 0:
                component_data.append(comp_rows)

        # If we have all components, sum them
        if len(component_data) == len(sw_components):
            # Get value columns (exclude 'usecol' and 'kper')
            value_cols = [col for col in component_data[0].columns if col not in ['usecol', 'kper']]

            # Sum the values
            summed_values = component_data[0][value_cols].copy()
            for comp_df in component_data[1:]:
                summed_values = summed_values.add(comp_df[value_cols].values, fill_value=0)

            # Create new rows for 'surface water'
            for orig_idx, (idx, row) in enumerate(summed_values.iterrows()):
                new_row = row.copy()
                new_row['usecol'] = 'surface water'
                new_row['kper'] = kper
                sw_data.append(new_row)

    # Append surface water data to budget_df
    if sw_data:
        sw_df = pd.DataFrame(sw_data)
        budget_df = pd.concat([budget_df, sw_df], ignore_index=True)

    unique_kpers = sorted(budget_df['kper'].unique())

    # Define desired order for usecols (exclude 'total', add 'surface water' after 'spring')
    desired_order = ['awanui', 'poukawa', 'spring', 'surface water', 'confined',
                     'recharge', 'stoss', 'inout', 'percentdiscrepancy']

    # Get all unique usecols from data
    all_usecols = set(budget_df['usecol'].unique())

    # Order usecols: first by desired order, then any remaining alphabetically
    unique_usecols = []
    for usecol in desired_order:
        if usecol in all_usecols:
            unique_usecols.append(usecol)
            all_usecols.remove(usecol)

    # Add any remaining usecols not in desired_order (excluding 'total')
    if all_usecols:
        remaining = sorted([u for u in all_usecols if u != 'total'])
        unique_usecols.extend(remaining)

    n_usecols = len(unique_usecols)

    def remove_outliers(data, method='iqr'):
        """Remove outliers using IQR method"""
        if len(data) == 0:
            return data
        q1 = np.percentile(data, 25)
        q3 = np.percentile(data, 75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        return data[(data >= lower_bound) & (data <= upper_bound)]

    # Organize subplots in rows with max 3 columns
    max_cols = 3
    n_cols = min(n_usecols, max_cols)
    n_rows = int(np.ceil(n_usecols / max_cols))

    # Create figure with subplots - share x-axis only
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 5 * n_rows),
                             sharex=True, sharey=False)

    # Flatten axes array for easier indexing
    if n_usecols == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if n_usecols > 1 else [axes]

    # Color scheme for kpers
    kper_colors = plt.cm.Set2(np.linspace(0, 1, len(unique_kpers)))

    # Track outlier statistics
    outlier_stats = []

    # Plot each usecol in its own subplot
    for idx, usecol in enumerate(unique_usecols):
        ax = axes[idx]

        # Get data for this usecol
        df_usecol = budget_df[budget_df['usecol'] == usecol]

        # Prepare data for boxplot - one boxplot per kper
        boxplot_data = []
        kper_labels = []

        for kper in unique_kpers:
            df_kper = df_usecol[df_usecol['kper'] == kper]
            value_cols = [col for col in df_kper.columns if col not in ['kper', 'usecol']]

            # Combine all values for this kper
            all_values = []
            for col in value_cols:
                values = df_kper[col].dropna().values
                all_values.extend(values)

            # Convert to array
            all_values = np.array(all_values)
            original_count = len(all_values)

            # Remove outliers
            if len(all_values) > 0:
                filtered_values = remove_outliers(all_values)
            else:
                filtered_values = all_values

            filtered_count = len(filtered_values)
            removed_count = original_count - filtered_count

            boxplot_data.append(filtered_values)
            kper_labels.append(f'kper {kper}')

            outlier_stats.append({
                'usecol': usecol,
                'kper': kper,
                'removed': removed_count,
                'original': original_count,
                'percent': removed_count / original_count * 100 if original_count > 0 else 0
            })

        # Create boxplot for this usecol
        bp = ax.boxplot(boxplot_data, labels=kper_labels, patch_artist=True)

        # Color the boxes by kper
        for patch, color in zip(bp['boxes'], kper_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Add individual points with jitter for each kper
        np.random.seed(42)  # For reproducible jitter
        for i, values in enumerate(boxplot_data):
            # Limit number of points to plot for clarity
            n_points = min(len(values), 200)
            if n_points < len(values):
                indices = np.random.choice(len(values), n_points, replace=False)
                values_to_plot = values[indices]
            else:
                values_to_plot = values

            if len(values_to_plot) > 0:
                x_pos = np.random.normal(i+1, 0.04, size=len(values_to_plot))
                ax.scatter(x_pos, values_to_plot, alpha=0.3, s=1, color='red')

        # Set labels and title
        ax.set_title(f'{usecol}', fontsize=11, fontweight='bold')
        ax.set_xlabel('Stress Period', fontsize=10)

        # Set y-axis label based on usecol type
        if usecol == 'percentdiscrepancy':
            ax.set_ylabel('Percent Discrepancy (%)', fontsize=10)
        else:
            ax.set_ylabel('Budget Value (m³/d)', fontsize=10)

        ax.grid(True, alpha=0.3, axis='y')

        # Rotate x-axis labels if needed
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    # Remove empty subplots if n_usecols doesn't fill all grid positions
    total_subplots = n_rows * n_cols
    if n_usecols < total_subplots:
        for idx in range(n_usecols, total_subplots):
            fig.delaxes(axes[idx])

    # Add main title
    fig.suptitle('Budget Components by Stress Period (Outliers Removed)',
                 fontsize=14, fontweight='bold', y=0.995)

    plt.tight_layout()

    # Save figure
    fig_fn = os.path.join(output, 'output.budget_boxplot_by_kper.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

    # Print summary of outlier removal
    print("\nOutlier removal summary (by kper):")
    print(f"{'Component':<30} {'Kper':<6} {'Removed':<10} {'Original':<10} {'Percent':<10}")
    print("-" * 70)
    for stat in outlier_stats:
        print(f"{stat['usecol']:<30} {stat['kper']:<6} {stat['removed']:<10} {stat['original']:<10} {stat['percent']:<10.1f}%")

def plot_parameter_histograms(data, subset_pe, output):
    """Plot parameter histograms by group - prior vs posterior."""
    par_data = data['par_data']
    pr_pe = data['pr_pe']
    pt_pe = subset_pe

    # Get parameter groups (excluding those starting with 's_1_' which are pilot points)
    all_groups = par_data['pargp'].unique()
    par_groups = [g for g in all_groups if not g.startswith('s_1_')]

    print(f"\nPlotting histograms for {len(par_groups)} parameter groups...")

    # Plot histograms for each parameter group
    n_groups = len(par_groups)
    n_cols = 3
    n_rows = int(np.ceil(n_groups / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    axes = axes.flatten() if n_groups > 1 else [axes]

    for idx, pargp in enumerate(sorted(par_groups)):
        ax = axes[idx]

        # Get parameters in this group
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()

        if len(pars_in_group) == 0:
            ax.text(0.5, 0.5, f'No parameters\nin {pargp}',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(pargp)
            continue

        # Filter to only parameters that exist in both prior and posterior
        pars_in_prior = [p for p in pars_in_group if p in pr_pe.columns]
        pars_in_posterior = [p for p in pars_in_group if p in pt_pe.columns]
        common_pars = list(set(pars_in_prior) & set(pars_in_posterior))

        if len(common_pars) == 0:
            ax.text(0.5, 0.5, f'No common parameters\nin {pargp}',
                   ha='center', va='center', transform=ax.transAxes)
            ax.set_title(pargp)
            continue

        # Get parameter transform
        transform = par_data.loc[common_pars[0]]['partrans']

        # Get values for all realizations and all parameters in this group
        prior_vals = pr_pe.loc[:, common_pars].values.flatten()
        posterior_vals = pt_pe.loc[:, common_pars].values.flatten()

        # Apply log transform if needed
        if transform == 'log':
            prior_vals = np.log10(prior_vals)
            posterior_vals = np.log10(posterior_vals)
            xlabel = f'log10(Parameter Value)'
        else:
            xlabel = 'Parameter Value'

        # Determine histogram bins based on combined data range
        all_vals = np.concatenate([prior_vals, posterior_vals])
        bins = np.linspace(np.percentile(all_vals, 1), np.percentile(all_vals, 99), 30)

        # Plot histograms
        ax.hist(prior_vals, bins=bins, alpha=0.5, color='grey',
               label='Prior', density=True, edgecolor='black', linewidth=0.5)
        ax.hist(posterior_vals, bins=bins, alpha=0.5, color='cornflowerblue',
               label='Posterior', density=True, edgecolor='black', linewidth=0.5)

        ax.set_xlabel(xlabel)
        ax.set_ylabel('Density')
        ax.set_title(f'{pargp}\n(n={len(common_pars)} params)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Remove empty subplots
    for idx in range(n_groups, len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    fig_fn = os.path.join(output, 'output.parameter_group_histograms.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

    # Print summary statistics
    print("\n" + "="*70)
    print("PARAMETER SUMMARY STATISTICS")
    print("="*70)
    for pargp in sorted(par_groups):
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
        if len(pars_in_group) == 0:
            continue

        # Filter to only parameters that exist in both prior and posterior
        pars_in_prior = [p for p in pars_in_group if p in pr_pe.columns]
        pars_in_posterior = [p for p in pars_in_group if p in pt_pe.columns]
        common_pars = list(set(pars_in_prior) & set(pars_in_posterior))

        if len(common_pars) == 0:
            continue

        prior_vals = pr_pe.loc[:, common_pars].values.flatten()
        posterior_vals = pt_pe.loc[:, common_pars].values.flatten()

        print(f"\n{pargp} ({len(common_pars)} parameters):")
        print(f"  Prior:     mean={np.mean(prior_vals):.4f}, std={np.std(prior_vals):.4f}")
        print(f"  Posterior: mean={np.mean(posterior_vals):.4f}, std={np.std(posterior_vals):.4f}")
        std_reduction = (1 - np.std(posterior_vals)/np.std(prior_vals))*100
        print(f"  Std reduction: {std_reduction:.1f}%")

def plot_1to1(data, output):
    """Plot 1:1 observed vs simulated."""
    pst = data['pst']

    fig = pst.plot(kind='1to1')
    fig_fn = os.path.join(output, 'output.1to1.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

def plot_1to1_boxplots(data, output):
    """Plot 1:1 observed vs simulated as boxplots showing ensemble uncertainty."""
    obs_data = data['obs_data']
    pr_oe = data['pr_oe']
    pt_oe = data['pt_oe']
    pst = data['pst']

    # Get observation groups
    obs_groups = obs_data['obgnme'].unique()

    # Filter to non-zero weighted observations
    nz_obs_names = obs_data[obs_data['weight'] > 0].index.tolist()

    if len(nz_obs_names) == 0:
        print("No non-zero weighted observations found, skipping boxplot 1:1")
        return

    # ObservationEnsemble is already a DataFrame-like object
    pr_oe_df = pr_oe
    pt_oe_df = pt_oe

    # Prepare data for plotting
    obs_plot = []
    prior_plot = []
    posterior_plot = []
    obs_names_plot = []

    for obs_name in nz_obs_names:
        try:
            obs_val = obs_data.loc[obs_name, 'obsval']

            # Get ensemble values using .loc
            prior_vals = pr_oe_df.loc[:, obs_name].values
            posterior_vals = pt_oe_df.loc[:, obs_name].values

            # Filter out NaN values
            prior_vals = prior_vals[~np.isnan(prior_vals)]
            posterior_vals = posterior_vals[~np.isnan(posterior_vals)]

            if len(prior_vals) > 0 and len(posterior_vals) > 0:
                obs_plot.append(obs_val)
                prior_plot.append(prior_vals)
                posterior_plot.append(posterior_vals)
                obs_names_plot.append(obs_name)
        except (KeyError, AttributeError):
            continue

    # Subsample if too many observations
    max_obs = 100
    if len(obs_plot) > max_obs:
        indices = np.random.choice(len(obs_plot), max_obs, replace=False)
        obs_plot = [obs_plot[i] for i in sorted(indices)]
        prior_plot = [prior_plot[i] for i in sorted(indices)]
        posterior_plot = [posterior_plot[i] for i in sorted(indices)]
        obs_names_plot = [obs_names_plot[i] for i in sorted(indices)]
        print(f"  Subsampled {max_obs} of {len(obs_names_plot)} observations for clarity")

    # PRIOR boxplot
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    # Create boxplot positions
    positions = np.arange(len(obs_plot))

    # Plot boxplots for prior
    bp = ax.boxplot(prior_plot, positions=positions, widths=0.6,
                     patch_artist=True, manage_ticks=False,
                     boxprops=dict(facecolor='lightgrey', alpha=0.7),
                     medianprops=dict(color='black', linewidth=1),
                     flierprops=dict(marker='o', markersize=2, alpha=0.3))

    # Plot observed values
    ax.scatter(positions, obs_plot, color='red', s=20, alpha=0.8,
              zorder=5, label='Observed', marker='D')

    # Add 1:1 line
    all_vals = np.concatenate([obs_plot] + prior_plot)
    min_val, max_val = np.percentile(all_vals[~np.isnan(all_vals)], [1, 99])
    ax.plot([0, len(obs_plot)-1], [min_val, max_val], 'k--',
           linewidth=1, alpha=0.5, label='1:1 line')

    ax.set_xlabel('Observation Index')
    ax.set_ylabel('Value')
    ax.set_title('Prior: Observed vs Simulated Ensemble')
    ax.legend(loc='best', fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.5, len(obs_plot)-0.5)

    # POSTERIOR boxplot
    ax = axes[1]

    # Plot boxplots for posterior
    bp = ax.boxplot(posterior_plot, positions=positions, widths=0.6,
                     patch_artist=True, manage_ticks=False,
                     boxprops=dict(facecolor='cornflowerblue', alpha=0.7),
                     medianprops=dict(color='black', linewidth=1),
                     flierprops=dict(marker='o', markersize=2, alpha=0.3))

    # Plot observed values
    ax.scatter(positions, obs_plot, color='red', s=20, alpha=0.8,
              zorder=5, label='Observed', marker='D')

    # Add 1:1 line
    all_vals = np.concatenate([obs_plot] + posterior_plot)
    min_val, max_val = np.percentile(all_vals[~np.isnan(all_vals)], [1, 99])
    ax.plot([0, len(obs_plot)-1], [min_val, max_val], 'k--',
           linewidth=1, alpha=0.5, label='1:1 line')

    ax.set_xlabel('Observation Index')
    ax.set_ylabel('Value')
    ax.set_title('Posterior: Observed vs Simulated Ensemble')
    ax.legend(loc='best', fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_xlim(-0.5, len(obs_plot)-0.5)

    plt.tight_layout()
    fig_fn = os.path.join(output, 'output.1to1_boxplots.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

    # Also create scatter-style 1:1 with uncertainty bands
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # Calculate percentiles for uncertainty bands
    prior_medians = [np.median(p) for p in prior_plot]
    prior_p25 = [np.percentile(p, 25) for p in prior_plot]
    prior_p75 = [np.percentile(p, 75) for p in prior_plot]

    posterior_medians = [np.median(p) for p in posterior_plot]
    posterior_p25 = [np.percentile(p, 25) for p in posterior_plot]
    posterior_p75 = [np.percentile(p, 75) for p in posterior_plot]

    # PRIOR scatter with error bars
    ax = axes[0]
    ax.errorbar(obs_plot, prior_medians,
                yerr=[np.array(prior_medians) - np.array(prior_p25),
                      np.array(prior_p75) - np.array(prior_medians)],
                fmt='o', color='grey', alpha=0.5, markersize=4,
                elinewidth=1, capsize=2, label='Prior (IQR)')

    # Add 1:1 line
    min_val = min(min(obs_plot), min(prior_medians))
    max_val = max(max(obs_plot), max(prior_medians))
    ax.plot([min_val, max_val], [min_val, max_val], 'k--',
           linewidth=2, label='1:1 line')

    ax.set_xlabel('Observed')
    ax.set_ylabel('Simulated (Median)')
    ax.set_title('Prior: 1:1 Plot with IQR')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    ax.axis('equal')

    # POSTERIOR scatter with error bars
    ax = axes[1]
    ax.errorbar(obs_plot, posterior_medians,
                yerr=[np.array(posterior_medians) - np.array(posterior_p25),
                      np.array(posterior_p75) - np.array(posterior_medians)],
                fmt='o', color='cornflowerblue', alpha=0.5, markersize=4,
                elinewidth=1, capsize=2, label='Posterior (IQR)')

    # Add 1:1 line
    min_val = min(min(obs_plot), min(posterior_medians))
    max_val = max(max(obs_plot), max(posterior_medians))
    ax.plot([min_val, max_val], [min_val, max_val], 'k--',
           linewidth=2, label='1:1 line')

    ax.set_xlabel('Observed')
    ax.set_ylabel('Simulated (Median)')
    ax.set_title('Posterior: 1:1 Plot with IQR')
    ax.legend(loc='best')
    ax.grid(alpha=0.3)
    ax.axis('equal')

    plt.tight_layout()
    fig_fn = os.path.join(output, 'output.1to1_scatter_uncertainty.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")


def extract_index(df, extract_dict):
    input_df = df.transpose()

    for col, func in extract_dict.items():
        input_df[col] = input_df.index.map(func).to_list()

    return input_df

def get_arr_list(df, col_names):
    # create a dummy array
    dummy_arr = np.full((47,57), np.nan)
    i_arrs = []
    for real in col_names:
        # fill array with values from df where i and j are indexes
        idf = df[['i', 'j', real]].copy()
        arr = dummy_arr.copy()
        for _, row in idf.iterrows():
            i = int(row['i']) - 1  # assuming i starts at 1
            j = int(row['j']) - 1  # assuming j starts at 1
            val = row[real]
            arr[i, j] = val
        # replace max_nan with np.nan
        max_nan = 1e+30
        arr = np.where(arr == max_nan, np.nan, arr)
        i_arrs.append(arr)
    return np.stack(i_arrs)


def get_arrs_stats(df, reals, is_kper=False):
    
    # get all the heads into an array which we can now call to plot
    arrs = {}
    if is_kper:
        key_iter = set(zip(df['kper'], df['kstp']))
    else:
        key_iter = set(zip([1], [1]))
    for i in key_iter:

        if is_kper:
            kper, kstp = i
            idf = df[(df['kper'] == kper) & (df['kstp'] == kstp)]
        else:
            idf = df

        i_stack = get_arr_list(idf, reals)
        arrs[i] = {
            'pr_stack': i_stack,
            'mean': np.mean(i_stack, axis=0),
            'std': np.std(i_stack, axis=0),
            'var': np.var(i_stack, axis=0),
            'min': np.min(i_stack, axis=0),
            'max': np.max(i_stack, axis=0),
        }
    return arrs

def plot_pname_array_statistics(data, output, pname_prefix, array_param=False):
    desired_cols = [col for col in data['param_headers'] if col.startswith(f'pname:{pname_prefix}')]

    if array_param:
        index_col=None
        extract_dict = {
            'i': extract_i,
            'j': extract_j,
        }
    else:
        index_col=0
        extract_dict = {
            'i': extract_idx1,
            'j': extract_idx2,
        }
    
    dfs = {
        'pe': pd.read_csv(data['pr_pe_fn'], usecols=desired_cols, index_col=index_col),
        'pt': pd.read_csv(data['pt_pe_fn'], usecols=desired_cols, index_col=index_col),
    }
    arrs = {}
    for name, rdf in dfs.items():
        rdf.reset_index(drop=True, inplace=True)
        reals = rdf.index.tolist()
        df = extract_index(rdf, extract_dict)
        arrs[name] = get_arrs_stats(df, reals)

        # plot & save
        plot_arr_stats(arrs[name], output, pname_prefix+f'_{name}')
    
    # percent change plot
    var_change = 100 * (1 - (arrs['pe'][(1,1)]['std'] / arrs['pt'][(1,1)]['std']))

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(var_change, cmap='bwr', vmin=-100, vmax=100)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title('Variance Percent Change (pt vs pe)')
    fig_fn = os.path.join(output, f'output.{pname_prefix}_percent_change.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()    


def plot_oname_array_statistics(data, subset_pe, output, oname_prefix):
    desired_cols = [col for col in data['headers'] if col.startswith(f'oname:{oname_prefix}')]
    
    # rdf = pd.read_csv(data[fn], usecols=desired_cols, index_col=0)
    
    extract_dict = {
        'i': extract_i,
        'j': extract_j,
        'kper': extract_kper,
        'kstp': extract_kstp,
    }

    pr_df = pd.read_csv(data['pr_oe_fn'], usecols=desired_cols, index_col=0)

    if subset_pe is not None:
        pe_df = subset_pe.loc[:, desired_cols]
    else:
        pe_df = pd.read_csv(data['pr_oe_fn'], usecols=desired_cols, index_col=0)

    dfs = {
        'pe': pr_df,
        'pt': pe_df,
    }
    
    arrs = {}
    for name, rdf in dfs.items():
        reals = rdf.index.tolist()

        df = extract_index(rdf, extract_dict)
    
        arrs[name] = get_arrs_stats(df, reals, is_kper=True)

        # plot & save
        plot_arr_stats(arrs[name], output, oname_prefix+f'_{name}')

    # percent change plot
    var_change = 100 * (1 - (arrs['pe'][(1,1)]['std'] / arrs['pt'][(1,1)]['std']))

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(var_change, cmap='bwr', vmin=-100, vmax=100)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title('Variance Percent Change (pt vs pe)')
    fig_fn = os.path.join(output, f'output.{oname_prefix}_var_percent_change.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()  


def plot_heads_oname_array_statistics(data, fn, output, oname_prefix, out_name):
    desired_cols = [col for col in data['headers'] if col.startswith(f'oname:{oname_prefix}')]
    
    rdf = pd.read_csv(data[fn], usecols=desired_cols, index_col=0)
    reals = rdf.index.tolist()
    extract_dict = {
        'i': extract_i,
        'j': extract_j,
        'kper': extract_hkper,
        'kstp': extract_hkstp,
    }

    df = extract_index(rdf, extract_dict)
    
    arrs = get_arrs_stats(df, reals, is_kper=True)

    # plot & save
    plot_arr_stats(arrs, output, out_name)



def plot_arr_stats(arrs, output, arr_name):

    # Determine available keys (handles both multiple kpers and single kper case)
    all_keys = list(arrs.keys())
    if isinstance(all_keys[0], tuple):
        # Filter for keys where kstp == 1 (second element of tuple)
        available_keys = sorted([key for key in all_keys if key[1] == 1])
    else:
        # Single key case (not a tuple)
        available_keys = sorted(all_keys)

    n_kpers = len(available_keys)

    # Create figure with shared axes and tight spacing
    # Adjust number of columns based on available kpers
    fig, axes = plt.subplots(nrows=4, ncols=n_kpers, figsize=(3.5*n_kpers, 12),
                             sharex=True, sharey=True,
                             gridspec_kw={'hspace': 0.05, 'wspace': 0.05})

    # Handle case where there's only 1 column (axes won't be 2D)
    if n_kpers == 1:
        axes = axes.reshape(-1, 1)

    # Calculate vmin/vmax for each row (stat type) across all available kpers
    mean_arrays = [arrs[key]['mean'] for key in available_keys]
    min_arrays = [arrs[key]['min'] for key in available_keys]
    max_arrays = [arrs[key]['max'] for key in available_keys]
    std_arrays = [arrs[key]['std'] for key in available_keys]

    # Calculate shared vmin/vmax for each row
    mean_vmin, mean_vmax = np.nanmin(mean_arrays), np.nanmax(mean_arrays)
    min_vmin, min_vmax = np.nanmin(min_arrays), np.nanmax(min_arrays)
    max_vmin, max_vmax = np.nanmin(max_arrays), np.nanmax(max_arrays)
    std_vmin, std_vmax = np.nanmin(std_arrays), np.nanmax(std_arrays)

    for idx, key in enumerate(available_keys):
        mean_arr = arrs[key]['mean']
        std_arr = arrs[key]['std']
        min_arr = arrs[key]['min']
        max_arr = arrs[key]['max']

        # Row 0: Mean
        im0 = axes[0, idx].imshow(mean_arr, cmap='viridis', vmin=mean_vmin, vmax=mean_vmax)
        cbar0 = plt.colorbar(im0, ax=axes[0, idx], fraction=0.046, pad=0.04)
        # Set title based on key (kper, kstp) or just single value
        if isinstance(key, tuple):
            axes[0, idx].set_title(f'kper {key[0]}', fontsize=10, pad=2)
        else:
            axes[0, idx].set_title(f'{key}', fontsize=10, pad=2)
        # Add text label showing mean value
        mean_val = np.nanmean(mean_arr)
        axes[0, idx].text(0.02, 0.98, f'μ={mean_val:.2e}',
                         transform=axes[0, idx].transAxes, fontsize=7,
                         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # Row 1: Min
        im1 = axes[1, idx].imshow(min_arr, cmap='plasma', vmin=min_vmin, vmax=min_vmax)
        cbar1 = plt.colorbar(im1, ax=axes[1, idx], fraction=0.046, pad=0.04)
        # Add text label showing min value
        min_val = np.nanmin(min_arr)
        axes[1, idx].text(0.02, 0.98, f'min={min_val:.2e}',
                         transform=axes[1, idx].transAxes, fontsize=7,
                         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # Row 2: Max
        im2 = axes[2, idx].imshow(max_arr, cmap='plasma', vmin=max_vmin, vmax=max_vmax)
        cbar2 = plt.colorbar(im2, ax=axes[2, idx], fraction=0.046, pad=0.04)
        # Add text label showing max value
        max_val = np.nanmax(max_arr)
        axes[2, idx].text(0.02, 0.98, f'max={max_val:.2e}',
                         transform=axes[2, idx].transAxes, fontsize=7,
                         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

        # Row 3: Std
        im3 = axes[3, idx].imshow(std_arr, cmap='magma', vmin=std_vmin, vmax=std_vmax)
        cbar3 = plt.colorbar(im3, ax=axes[3, idx], fraction=0.046, pad=0.04)
        # Add text label showing mean std value
        std_val = np.nanmean(std_arr)
        axes[3, idx].text(0.02, 0.98, f'σ={std_val:.2e}',
                         transform=axes[3, idx].transAxes, fontsize=7,
                         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    # Add row labels on the left
    axes[0, 0].set_ylabel('Mean\n\ni', fontsize=9)
    axes[1, 0].set_ylabel('Min\n\ni', fontsize=9)
    axes[2, 0].set_ylabel('Max\n\ni', fontsize=9)
    axes[3, 0].set_ylabel('Std\n\ni', fontsize=9)

    # Add x-labels only to bottom row
    for idx in range(n_kpers):
        axes[3, idx].set_xlabel('j', fontsize=9)

    # Remove tick labels except for outer edges
    for i in range(4):
        for j in range(n_kpers):
            if i < 3:  # Not bottom row
                axes[i, j].set_xticklabels([])
            if j > 0:  # Not left column
                axes[i, j].set_yticklabels([])

    fig_fn = os.path.join(output, f'arrs.pt_{arr_name}.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")


def plot_pdc_budget_confined(data, output_dir):
    """Plot budget term histograms for all stress periods in a single figure."""
    desired_cols = [col for col in data['headers'] if (col.startswith(f'oname:budget')) & (col.split('_')[2].split(':')[1] == 'confined')]

    if len(desired_cols) == 0:
        print("No confined budget observations found, skipping plot")
        return

    rdf = pd.read_csv(data['pr_oe_fn'], usecols=desired_cols)
    noise = pd.read_csv(data['noise_fn'], usecols=desired_cols)

    # Sort columns by kper for consistent ordering
    desired_cols = sorted(desired_cols)
    n_cols = len(desired_cols)

    # Create figure with subplots (2 rows x 2 cols for 4 stress periods)
    n_rows = int(np.ceil(n_cols / 2))
    fig, axes = plt.subplots(n_rows, 2, figsize=(12, 4*n_rows))
    axes = axes.flatten() if n_cols > 1 else [axes]

    # Plot histogram for each budget term
    for idx, col in enumerate(desired_cols):
        ax = axes[idx]

        # Prior budget values
        prior_vals = rdf[col].values
        prior_vals = prior_vals[~np.isnan(prior_vals)]

        # Noise budget values
        noise_vals = noise[col].values
        noise_vals = noise_vals[~np.isnan(noise_vals)]

        # Determine bins
        all_vals = np.concatenate([prior_vals, noise_vals])
        bins = np.linspace(np.percentile(all_vals, 1), np.percentile(all_vals, 99), 30)

        # Plot prior histogram
        ax.hist(prior_vals, bins=bins, alpha=0.5, color='grey',
                label='Prior', density=True, edgecolor='black', linewidth=0.5)

        # Plot noise histogram
        ax.hist(noise_vals, bins=bins, alpha=0.5, color='orange',
                label='Obs+Noise', density=True, edgecolor='darkorange', linewidth=0.5)

        ax.set_xlabel('Budget Term Value', fontsize=10)
        ax.set_ylabel('Density', fontsize=10)
        ax.set_title(f'kper: {col[-1]}', fontsize=11, fontweight='bold')
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3)

        # Add summary statistics as text
        prior_mean = np.mean(prior_vals)
        noise_mean = np.mean(noise_vals)

        # add lines to plot of mean
        ax.axvline(prior_mean, color='black', linestyle='--', linewidth=1.5, label='Prior Mean')
        ax.axvline(noise_mean, color='darkorange', linestyle='--', linewidth=1.5, label='Obs+Noise Mean')

        ax.text(0.98, 0.98, f'Prior μ={prior_mean:.2e}\nNoise μ={noise_mean:.2e}',
                transform=ax.transAxes, fontsize=8,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    # Remove empty subplots if n_cols is odd
    if n_cols < len(axes):
        for idx in range(n_cols, len(axes)):
            fig.delaxes(axes[idx])

    plt.tight_layout()
    fig_fn = os.path.join(output_dir, 'output.pdc_budget_confined.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")


def plot_histograms_spring_flux(data, subset_oe, output, col_startswith='ts-flux'):
    """Plot time series observations - prior vs posterior vs truth."""
    headers = data['headers']
    pr_oe_fn = data['pr_oe_fn']
    # pt_oe_fn = data['pt_oe_fn']

    # Get time series head observations
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("No time series head observations found, skipping plot")
        return

    pre_oe_ts = pd.read_csv(pr_oe_fn, usecols=desired_cols)
    pt_oe_ts = subset_oe.loc[:, desired_cols]

    ts_df = format_tsdf(pre_oe_ts)
    ts_pt_df = format_tsdf(pt_oe_ts)

    usecol = ts_df['usecol'].unique()[0]

    # Prep data
    df_sub = ts_df[(ts_df['usecol'] == usecol)].sort_values('time').replace(-999, np.nan)
    df_pt_sub = ts_pt_df[(ts_pt_df['usecol'] == usecol)].sort_values('time').replace(-999, np.nan)

    value_cols = [col for col in df_sub.columns if col not in ['time', 'kper', 'kstp', 'usecol']]
    # Ensure columns exist in posterior and noise dataframes
    value_cols_pt = [col for col in df_pt_sub.columns if col not in ['time', 'kper', 'kstp', 'usecol']]

    fig, axes = plt.subplots(nrows=1, ncols= 4, figsize=(6, 2))

    kpers = [1, 2, 3, 4]

    for i, kper in enumerate(kpers):
        kper_mask = df_sub['kper'] == kper
        plt_pr = df_sub[kper_mask]
        plt_pt = df_pt_sub[kper_mask]

        # Determine shared bins for consistent plotting across both histograms
        pr_values = plt_pr[value_cols].values.flatten()
        pt_values = plt_pt[value_cols_pt].values.flatten()
        # all_vals = np.concatenate([pr_values, pt_values])

        if len(pt_values) > 0:
            min = np.percentile(pt_values, 5)
            max = np.percentile(pt_values, 95)
            bins = np.linspace(min, max, 20)
        else:
            min = 0
            max = 1
            bins = 20
        #     
        # else:
        #     bins = 15

        axes[i].hist(
            pr_values, bins=20, alpha=0.5, color='lightgrey', 
            density=True, edgecolor='grey', linewidth=0.8,
            label='Prior')
        axes[i].hist(
            pt_values, bins=bins, alpha=0.5, color='blue',
            density=True, edgecolor='darkblue', linewidth=0.8,
            label='Posterior')

        # Add mean line for "present" (kper 1-2)
        mean_present = np.mean(pt_values)
        axes[i].axvline(mean_present, color='purple', linestyle='-', linewidth=2, label='Mean')

        axes[i].set_xlabel(f'Flux (m\u00b3/d)', fontsize=8)
        # ax_hist.set_title('kper 1-2', fontsize=9)
        axes[i].set_xlim(min, max)
        axes[i].tick_params(axis='both', labelsize=8)
        axes[i].grid(alpha=0.3)
        plt.setp(axes[i].get_yticklabels(), visible=False)

        # Add legend only to first subplot
        if i == 0:
            axes[i].legend(fontsize=6, loc='upper left')

    plt.tight_layout()
    fig_fn = os.path.join(output, 'prediction.spring_flux_kper.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

def plot_timeseries_spring_flux(data, subset_oe, output, col_startswith='ts-flux'):
    """Plot time series observations - prior vs posterior vs truth."""
    headers = data['headers']
    pr_oe_fn = data['pr_oe_fn']

    # Get time series head observations
    desired_cols = [col for col in headers if col.startswith(f'oname:{col_startswith}')]

    if len(desired_cols) == 0:
        print("No time series head observations found, skipping plot")
        return

    pre_oe_ts = pd.read_csv(pr_oe_fn, usecols=desired_cols)
    pt_oe_ts = subset_oe.loc[:, desired_cols]

    ts_df = format_tsdf(pre_oe_ts)
    ts_pt_df = format_tsdf(pt_oe_ts)

    usecol = ts_df['usecol'].unique()[0]
    n = 2 # two predicted periods

    # Prep data
    df_sub = ts_df[(ts_df['usecol'] == usecol)].sort_values('time').replace(-999, np.nan)
    df_pt_sub = ts_pt_df[(ts_pt_df['usecol'] == usecol)].sort_values('time').replace(-999, np.nan)

    value_cols = [col for col in df_sub.columns if col not in ['time', 'kper', 'kstp', 'usecol']]
    # Ensure columns exist in posterior and noise dataframes
    value_cols_pt = [col for col in df_pt_sub.columns if col not in ['time', 'kper', 'kstp', 'usecol']]

    fig, ax = plt.subplots(nrows=n, ncols=1, figsize=(6, 2*n), sharex=True)

    kpers = [2, 4]

    for i, kper in enumerate(kpers):
        kper_mask = df_sub['kper'] == kper
        plt_pr = df_sub[kper_mask]
        plt_pt = df_pt_sub[kper_mask]

        # Plot ensemble members
        for col in value_cols:
            y = plt_pr[col].copy() if col in plt_pr.columns else None

            if y is None and y_pt is None:
                continue

            # Mask outliers
            if y is not None:
                ax[i].plot(np.arange(0,len(y)), y, alpha=0.3, color='grey', linewidth=0.5)

    for i, kper in enumerate(kpers):
        kper_mask = df_sub['kper'] == kper
        plt_pt = df_pt_sub[kper_mask]

        # Plot ensemble members
        for col in value_cols:
            y_pt = plt_pt[col].copy() if col in plt_pt.columns else None

            if y is None and y_pt is None:
                continue

            # Mask outliers
            if y is not None:
                ax[i].plot(np.arange(0,len(y)), y, alpha=0.3, color='grey', linewidth=0.5)

            if y_pt is not None:
                ax[i].plot(np.arange(0,len(y_pt)), y_pt, alpha=0.3, color='blue', linewidth=0.8)

        if i == 0:
            ax[i].set_title(f'Present', loc='left')
        else:
            ax[i].set_title(f'Pre-European', loc='left')
        ax[i].set_ylabel('Spring flux (m3/d)')
        ax[i].grid(alpha=0.3)

        if i == 0:
            ax[i].legend(['Prior', 'Posterior'], fontsize=8)

        # Only show x-label on bottom plot
        if i == n - 1:
            ax[i].set_xlabel('Days')
        else:
            ax[i].set_xlabel('')

    # present day
    ax[0].set_ylim(np.min(df_pt_sub[value_cols_pt].values), np.max(df_pt_sub[value_cols_pt].values))

    plt.tight_layout()
    fig_fn = os.path.join(output, 'prediction.timeseries_spring_flux.png')
    plt.savefig(fig_fn, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {fig_fn}")

def plot_parameter_distributions(data, output_dir, subset_oe=None):
    
    par_data = data['par_data']
    pr_pe = data['pr_pe']
    pt_pe = data['pt_pe']

    # Get parameter groups (excluding those starting with 's_1_')
    all_groups = par_data['pargp'].unique()
    par_groups = [g for g in all_groups if not g.startswith('s_1_')]
    print(f"\nParameter groups to plot ({len(par_groups)}):")
    for g in sorted(par_groups):
        print(f"  - {g}")

    # Filter pt_reals to only include valid realization indices
    if subset_oe:
        valid_pt_reals = [r for r in pt_reals_list if r < len(pt_pe)]
        if len(valid_pt_reals) < len(pt_reals_list):
            print(f"\nWarning: Only {len(valid_pt_reals)} of {len(pt_reals_list)} requested realizations exist")
            print(f"Valid realizations: {valid_pt_reals}")

        # Create subset of posterior for highlighted realizations
        pt_subset = pt_pe.iloc[valid_pt_reals]

    # Plot histograms for each parameter group
    n_groups = len(par_groups)
    n_cols = 3
    n_rows = int(np.ceil(n_groups / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    axes = axes.flatten() if n_groups > 1 else [axes]

    for idx, pargp in enumerate(sorted(par_groups)):
        ax = axes[idx]

        # Get parameters in this group
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
        # get parameter transform
        transform = par_data.loc[pars_in_group]['partrans'][0]

        if len(pars_in_group) == 0:
            ax.text(0.5, 0.5, f'No parameters\nin {pargp}',
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title(pargp)
            continue

        # Get values for all realizations and all parameters in this group
        prior_vals = pr_pe.loc[:, pars_in_group].values.flatten()
        posterior_vals = pt_pe.loc[:, pars_in_group].values.flatten()

        if subset:
            subset_vals = pt_subset.loc[:, pars_in_group].values.flatten()

        if transform == 'log':
            prior_vals = np.log(prior_vals)
            posterior_vals = np.log(posterior_vals)

            if subset:
                subset_vals = np.log(subset_vals)

        # Determine histogram bins based on combined data range
        all_vals = np.concatenate([prior_vals, posterior_vals])
        bins = np.linspace(np.percentile(all_vals, 1), np.percentile(all_vals, 99), 30)

        # Plot histograms
        ax.hist(prior_vals, bins=bins, alpha=0.5, color='grey',
                label='Prior', density=True, edgecolor='black', linewidth=0.5)
        ax.hist(posterior_vals, bins=bins, alpha=0.5, color='cornflowerblue',
                label='Posterior', density=True, edgecolor='black', linewidth=0.5)
        if subset:
            ax.hist(subset_vals, bins=bins, alpha=0.7, color='darkblue',
                    label=f'Subset (n={len(valid_pt_reals)})', density=True,
                    edgecolor='black', linewidth=0.5)

        ax.set_xlabel('Parameter Value')
        ax.set_ylabel('Density')
        ax.set_title(f'{pargp}\n(n={len(pars_in_group)} params)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Remove empty subplots
    for idx in range(n_groups, len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'output.parameter_group_histograms.png'),
                dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {os.path.join(output_dir, 'output.parameter_group_histograms.png')}")

    # Optional: Print summary statistics
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    for pargp in sorted(par_groups):
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
        if len(pars_in_group) == 0:
            continue

        prior_vals = pr_pe.loc[:, pars_in_group].values.flatten()
        posterior_vals = pt_pe.loc[:, pars_in_group].values.flatten()

        print(f"\n{pargp} ({len(pars_in_group)} parameters):")
        print(f"  Prior:     mean={np.mean(prior_vals):.4f}, std={np.std(prior_vals):.4f}")
        print(f"  Posterior: mean={np.mean(posterior_vals):.4f}, std={np.std(posterior_vals):.4f}")
        print(f"  Std reduction: {(1 - np.std(posterior_vals)/np.std(prior_vals))*100:.1f}%")


def plot_parameter_distributions_subset(data, output_dir, subset_pe=None):
    
    par_data = data['par_data']
    pr_pe = data['pr_pe']
    pt_pe = subset_pe

    # Get parameter groups (excluding those starting with 's_1_')
    all_groups = par_data['pargp'].unique()
    par_groups = [g for g in all_groups if not g.startswith('s_1_')]
    print(f"\nParameter groups to plot ({len(par_groups)}):")
    for g in sorted(par_groups):
        print(f"  - {g}")

    # Plot histograms for each parameter group
    n_groups = len(par_groups)
    n_cols = 3
    n_rows = int(np.ceil(n_groups / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 4*n_rows))
    axes = axes.flatten() if n_groups > 1 else [axes]

    for idx, pargp in enumerate(sorted(par_groups)):
        ax = axes[idx]

        # Get parameters in this group
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
        # get parameter transform
        transform = par_data.loc[pars_in_group]['partrans'][0]

        if len(pars_in_group) == 0:
            ax.text(0.5, 0.5, f'No parameters\nin {pargp}',
                    ha='center', va='center', transform=ax.transAxes)
            ax.set_title(pargp)
            continue

        # Get values for all realizations and all parameters in this group
        prior_vals = pr_pe.loc[:, pars_in_group].values.flatten()
        posterior_vals = pt_pe.loc[:, pars_in_group].values.flatten()

        if transform == 'log':
            prior_vals = np.log(prior_vals)
            posterior_vals = np.log(posterior_vals)

        # Determine histogram bins based on combined data range
        all_vals = np.concatenate([prior_vals, posterior_vals])
        bins = np.linspace(np.percentile(all_vals, 1), np.percentile(all_vals, 99), 30)

        # Plot histograms
        ax.hist(prior_vals, bins=bins, alpha=0.5, color='grey',
                label='Prior', density=True, edgecolor='black', linewidth=0.5)
        ax.hist(posterior_vals, bins=bins, alpha=0.5, color='cornflowerblue',
                label='Posterior', density=True, edgecolor='black', linewidth=0.5)

        ax.set_xlabel('Parameter Value')
        ax.set_ylabel('Density')
        ax.set_title(f'{pargp}\n(n={len(pars_in_group)} params)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)

    # Remove empty subplots
    for idx in range(n_groups, len(axes)):
        fig.delaxes(axes[idx])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'output.parameter_group_histograms.png'),
                dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {os.path.join(output_dir, 'output.parameter_group_histograms.png')}")

    # Optional: Print summary statistics
    print("\n" + "="*70)
    print("SUMMARY STATISTICS")
    print("="*70)
    for pargp in sorted(par_groups):
        pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
        if len(pars_in_group) == 0:
            continue

        prior_vals = pr_pe.loc[:, pars_in_group].values.flatten()
        posterior_vals = pt_pe.loc[:, pars_in_group].values.flatten()

        print(f"\n{pargp} ({len(pars_in_group)} parameters):")
        print(f"  Prior:     mean={np.mean(prior_vals):.4f}, std={np.std(prior_vals):.4f}")
        print(f"  Posterior: mean={np.mean(posterior_vals):.4f}, std={np.std(posterior_vals):.4f}")
        print(f"  Std reduction: {(1 - np.std(posterior_vals)/np.std(prior_vals))*100:.1f}%")

