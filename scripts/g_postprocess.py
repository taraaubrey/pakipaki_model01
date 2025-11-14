"""
Post-processing script to generate figures from PEST++ IES results.
Adapted from notebooks/figure.ipynb

Usage:
    # For history matching (posterior) results:
    python scripts/g_postprocess_figures.py

    # For prior-only results:
    python scripts/g_postprocess_figures.py --prior-only
"""

import os
from g_helpers_run4 import *


def main():
    """Main post-processing function."""

    run_name = 'local_run30'
    last_iter = 9

    # Setup paths
    paths = setup_paths(run_name, m_d_str='master_ies', m_d_prior='master_ies')

    check_dir = paths['m_d']
    if not os.path.exists(check_dir):
        print(f"ERROR: Master IES directory not found: {check_dir}")
        print("Please run PEST++ IES first (d_pest_IES_HM.py)")
        return

    # Load data
    print("Loading data...")
    data = load_data(paths, run_name, last_iter=last_iter, use_prior_only=False, restart=False)
    
    output_dir = paths['output']

    # subset_indexes, means = subset_realizations_by_mean_heads(data, 'arr-h')
    subset_indexes, phis = subset_realizations_by_phi(data, min_threshold=300)
    # subset_indexes, means = subset_realizations_by_mean_heads(data, 'arr-h', from_iter=9)
    
    subset_oe_fn = os.path.join(paths['output'], 'subset_oe.csv')
    subset_pe_fn = os.path.join(paths['output'], 'subset_pe.csv')
    if not os.path.exists(subset_oe_fn):
        subset_oe, subset_pe = get_subset_ensembles(paths, data, run_name, subset_indexes, head_filter=True)
    else:
        subset_oe = pd.read_csv(subset_oe_fn)
        subset_pe = pd.read_csv(subset_pe_fn)

    # plot_oname_array_statistics(data, subset_oe, output_dir, 'arr-h')

    # Generate figures
    print("\n--- Generating Figures ---\n")

    # These plots only make sense for history matching runs
    print("1. Phi comparison plots...")
    plot_phi_comparison(data, output_dir)
    plot_phi_boxplot(data, output_dir)

    # plot parameter distributions
    # plot_parameter_distributions(data, output_dir, subset_pe)
    plot_parameter_distributions_subset(data, output_dir, subset_pe)
    
    print("\n2. Prediction plots...")
    plot_histograms_spring_flux(data, subset_oe, output_dir, col_startswith='ts-sprflux')
    plot_timeseries_spring_flux(data, subset_oe, output_dir, col_startswith='ts-sprflux')
    
    print("\n3. Observation plots...")
    plot_timeseries_obs(
        data, subset_oe, output_dir, 
        truth_file='truth/output.sample_heads.truth.csv', 
        col_startswith='ts-heads')
    
    plot_budgets(data, subset_oe, output_dir, col_startswith='budget')
    plot_budget_histogram(data, output_dir, col_startswith='budget')

    # print("\n3. Parameter histogram plots...")
    plot_parameter_histograms(data, subset_pe, output_dir)

    print("\n4. Parameter array statistics (posterior)...")
    # plot_pname_array_statistics(data, output_dir, 'ghbconf-cond-gr')
    # plot_pname_array_statistics(data, output_dir, 'ghbconf-head-hg')
    # plot_pname_array_statistics(data, output_dir, 'ghbaw-cond-gr')
    # plot_pname_array_statistics(data, output_dir, 'ghbaw-head-hg')
    # plot_pname_array_statistics(data, output_dir, 'npfklayer1-gr', array_param=True)

    print("\n5. Observation array statistics (posterior)...")
    # plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'confq', 'prior_confq')
    plot_oname_array_statistics(data, subset_oe, output_dir, 'arr-h')
    plot_oname_array_statistics(data, output_dir, 'arr-confq')
    plot_oname_array_statistics(data, 'pt_oe_fn', output_dir, 'awq', 'posterior_awq')

    print("\n6. 1:1 observed vs simulated plot...")
    plot_1to1(data, output_dir)

    print("\n7. 1:1 boxplot and scatter with uncertainty...")
    plot_1to1_boxplots(data, output_dir)





    # Prior-only mode: only generate array statistics
    print("1. Parameter array statistics (prior)...")
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'ghb-conf-cond-gr', 'ghb-conf-cond-gr')
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'ghb-conf-head-gr', 'ghb-conf-head-gr')
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'npfklayer1-fngr', 'npfklayer1-fngr', array_param=True)
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'npfklayer1-pp', 'npfklayer1-pp', array_param=True)
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'rch-fngr', 'rch-fngr', array_param=True)
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'rch-pp', 'rch-pp', array_param=True)
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'stosslayer1-pp', 'stosslayer1-pp', array_param=True)
    plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'stosslayer1-fngr', 'stosslayer1-fngr', array_param=True)

    print("\n2. Observation array statistics (prior)...")
    plot_heads_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'h-lyr', 'prior_heads')
    plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'confq', 'prior_confq')
    plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'awq', 'prior_awq')

    print("\n3. Check for prior data conflicts...")
    plot_pdc_budget_confined(data, output_dir)
        

    print(f"\n{'='*70}")
    print(f"POST-PROCESSING COMPLETE!")
    print(f"Figures saved to: {output_dir}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
