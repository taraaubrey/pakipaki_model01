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

def plot_obs_hist(pr_oe, pt_oe):
    for col in pr_oe.columns:
        # plot the prior and posterior ensemble distributions
        fig, ax = plt.subplots(1, 1)
        pr_oe.loc[:, col].hist(ax=ax, fc="0.5", ec="none", alpha=0.5, density=True)  # prior grey
        pt_oe.loc[:, col].hist(ax=ax, fc="b", ec="none", alpha=0.5, density=True)  # posterior blue
        ax.set_title(col)
        plt.show()

        return