"""
Plot parameter group histograms comparing prior, posterior, and subset of realizations.
Excludes parameter groups starting with 's_1_'.

Usage in Jupyter:
    %run scripts/plot_parameter_histograms.py
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyemu

# Configuration - UPDATE THESE PATHS
run_name = 'local_run1'  # Change to your run name
m_d = f'models/{run_name}/pest/master_ies'
m_d_pr = f'models/{run_name}/pest/master_ies_prior'
output_dir = f'models/{run_name}/figures'

# Define which posterior realizations to highlight
# pt_reals_list = [0, 1, 2, 3, 4]  # Example: highlight first 5 realizations
pt_reals_list = list(range(10))  # Example: highlight first 10 realizations

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Load PST file
pst = pyemu.Pst(os.path.join(m_d, f"{run_name}.pst"))

# Load parameter data
par_data = pst.parameter_data

# Determine the final iteration number
par_files = [f for f in os.listdir(m_d) if f.endswith('.par.csv') and not 'rejected' in f]
iterations = sorted([int(f.split('.')[1]) for f in par_files if f.split('.')[1].isdigit()])
final_iter = max(iterations) if iterations else 0

print(f"Loading parameter ensembles...")
print(f"Prior from: {m_d_pr}")
print(f"Posterior (iteration {final_iter}) from: {m_d}")

# Load parameter ensembles
try:
    prior_pe = pyemu.ParameterEnsemble.from_csv(
        pst=pst,
        filename=os.path.join(m_d_pr, f"{run_name}.0.par.csv")
    )
    print(f"Prior ensemble loaded: {prior_pe.shape}")
except Exception as e:
    print(f"Error loading prior ensemble: {e}")
    raise

try:
    posterior_pe = pyemu.ParameterEnsemble.from_csv(
        pst=pst,
        filename=os.path.join(m_d, f"{run_name}.{final_iter}.par.csv")
    )
    print(f"Posterior ensemble loaded: {posterior_pe.shape}")
except Exception as e:
    print(f"Error loading posterior ensemble: {e}")
    raise

# Get parameter groups (excluding those starting with 's_1_')
all_groups = par_data['pargp'].unique()
par_groups = [g for g in all_groups if not g.startswith('s_1_')]
print(f"\nParameter groups to plot ({len(par_groups)}):")
for g in sorted(par_groups):
    print(f"  - {g}")

# Filter pt_reals to only include valid realization indices
valid_pt_reals = [r for r in pt_reals_list if r < len(posterior_pe)]
if len(valid_pt_reals) < len(pt_reals_list):
    print(f"\nWarning: Only {len(valid_pt_reals)} of {len(pt_reals_list)} requested realizations exist")
    print(f"Valid realizations: {valid_pt_reals}")

# Create subset of posterior for highlighted realizations
pt_subset = posterior_pe.iloc[valid_pt_reals]

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

    # Get values for all realizations and all parameters in this group
    prior_vals = prior_pe.loc[:, pars_in_group].values.flatten()
    posterior_vals = posterior_pe.loc[:, pars_in_group].values.flatten()
    subset_vals = pt_subset.loc[:, pars_in_group].values.flatten()

    # Determine histogram bins based on combined data range
    all_vals = np.concatenate([prior_vals, posterior_vals])
    bins = np.linspace(np.percentile(all_vals, 1), np.percentile(all_vals, 99), 30)

    # Plot histograms
    ax.hist(prior_vals, bins=bins, alpha=0.5, color='grey',
            label='Prior', density=True, edgecolor='black', linewidth=0.5)
    ax.hist(posterior_vals, bins=bins, alpha=0.5, color='cornflowerblue',
            label='Posterior', density=True, edgecolor='black', linewidth=0.5)
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
plt.savefig(os.path.join(output_dir, 'parameter_group_histograms.png'),
            dpi=300, bbox_inches='tight')
print(f"\nFigure saved to: {os.path.join(output_dir, 'parameter_group_histograms.png')}")
plt.show()

# Optional: Print summary statistics
print("\n" + "="*70)
print("SUMMARY STATISTICS")
print("="*70)
for pargp in sorted(par_groups):
    pars_in_group = par_data[par_data['pargp'] == pargp].index.tolist()
    if len(pars_in_group) == 0:
        continue

    prior_vals = prior_pe.loc[:, pars_in_group].values.flatten()
    posterior_vals = posterior_pe.loc[:, pars_in_group].values.flatten()

    print(f"\n{pargp} ({len(pars_in_group)} parameters):")
    print(f"  Prior:     mean={np.mean(prior_vals):.4f}, std={np.std(prior_vals):.4f}")
    print(f"  Posterior: mean={np.mean(posterior_vals):.4f}, std={np.std(posterior_vals):.4f}")
    print(f"  Std reduction: {(1 - np.std(posterior_vals)/np.std(prior_vals))*100:.1f}%")
