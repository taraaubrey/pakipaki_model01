import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_phi_comparison(f_act,f_meas, output):
    fig, axes = plt.subplots(1, 2, sharey=True, figsize=(10,3.5))
    # left
    ax = axes[0]
    phi = pd.read_csv(f_act, index_col=0)
    phi.index = phi.total_runs
    phi.iloc[:,6:].apply(np.log10).plot(legend=False,lw=0.5,color='k', ax=ax)
    ax.set_title(r'Actual $\Phi$')
    ax.set_ylabel(r'log $\Phi$')
    # right
    ax = axes[-1]
    phi = pd.read_csv(f_meas,index_col=0)
    phi.index = phi.total_runs
    phi.iloc[:,6:].apply(np.log10).plot(legend=False,lw=0.2,color='r', ax=ax)
    ax.set_title(r'Measured+Noise $\Phi$')
    fig.tight_layout()

    # save
    fig_fn = os.path.join(output, 'output.phi_comparison.png')
    plt.savefig(fig_fn, dpi=300)

    plt.show()

    return

def plot_phi_hist(pr_oe, pt_oe):
    fig,ax = plt.subplots(1,1)
    pr_oe.phi_vector.apply(np.log10).hist(ax=ax,fc="0.5",ec="none",alpha=0.5,density=False) # prior grey
    pt_oe.phi_vector.apply(np.log10).hist(ax=ax,fc="b",ec="none",alpha=0.5,density=False) # posterior blue
    _ = ax.set_xlabel("$log_{10}\\phi$")

    return

def plot_obs_budget_hist(pr_oe, pt_oe, headers, obs_data, output):
    grp_start='oname:budget'
    # Decide which columns you want to load
    desired_cols1 = [col for col in headers if (col.startswith(grp_start) and col.endswith('1'))]
    desired_cols3 = [col for col in headers if (col.startswith(grp_start) and col.endswith('3'))]

    for desired_cols in [desired_cols1, desired_cols3]:
        kper = int(desired_cols[0].split('_')[-1].split(':')[-1])  # extract kper from the first column name
        # Calculate subplot dimensions
        n_cols = 4
        n_plots = len(desired_cols)
        n_rows = (n_plots + n_cols - 1) // n_cols  # ceiling division
        
        # Create subplots
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4 * n_rows))
        
        # Handle case where we have only one row
        if n_rows == 1:
            axes = axes.reshape(1, -1)
        elif n_plots == 1:
            axes = np.array([[axes]])
        
        # Flatten axes for easy iteration
        axes_flat = axes.flatten()
        
        for i, col in enumerate(desired_cols):
            ax = axes_flat[i]
            ax_b = ax.twinx()
            
            # plot the prior and posterior ensemble distributions
            pr_oe.loc[:, col].hist(ax=ax, fc="0.5", ec="none", alpha=0.5, density=True)  # prior grey
            pt_oe.loc[:, col].hist(ax=ax_b, fc="b", ec="none", alpha=0.5, density=True)  # posterior blue
            # plot x line
            ax.axvline(obs_data.loc[col, 'obsval'], color='red')
            ax.set_title(col, fontsize=10)
        
        # Hide empty subplots
        for i in range(n_plots, len(axes_flat)):
            axes_flat[i].set_visible(False)
        
        plt.tight_layout()

        # save
        fig_fn = os.path.join(output, f'output.budget_hist_kper{kper}.png')
        plt.savefig(fig_fn, dpi=300)

    return

def plot_obs_tsheads_linesplot(pr_oe, pt_oe, headers, obs_data, output):
    grp_start='oname:ts-heads'
    # Decide which columns you want to load
    desired_cols1 = [col for col in headers if (col.startswith(grp_start))]
    desired_cols3 = [col for col in headers if (col.startswith(grp_start) and col.endswith('3'))]