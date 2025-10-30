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
import argparse

from g_helpers import *
from setup import MODEL_NAME


def main():
    """Main post-processing function."""
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Generate figures from PEST++ results')
    parser.add_argument('--prior-only', action='store_true',
                       help='Use prior data only from master_ies_prior directory')
    args = parser.parse_args()

    # args.prior_only = True
    # Use MODEL_NAME from setup.py
    run_name = MODEL_NAME
    mode_str = "PRIOR-ONLY" if args.prior_only else "HISTORY MATCHING"
    print(f"\n{'='*70}")
    print(f"POST-PROCESSING FIGURES FOR: {run_name} ({mode_str})")
    print(f"{'='*70}\n")

    # Setup paths
    paths = setup_paths(run_name)

    # Check if appropriate directory exists
    if args.prior_only:
        check_dir = paths['m_d_pr']
        if not os.path.exists(check_dir):
            print(f"ERROR: Master IES Prior directory not found: {check_dir}")
            print("Please run PEST++ IES Prior first (c_pest_IES_prior.py)")
            return
    else:
        check_dir = paths['m_d']
        if not os.path.exists(check_dir):
            print(f"ERROR: Master IES directory not found: {check_dir}")
            print("Please run PEST++ IES first (d_pest_IES_HM.py)")
            return

    # Load data
    print("Loading data...")
    data = load_data(paths, run_name, use_prior_only=args.prior_only)
    output_dir = paths['output']

    # Add suffix to output filenames for prior-only mode
    if args.prior_only:
        print("\n*** PRIOR-ONLY MODE: Only observation and parameter array statistics will be generated ***\n")

    # Generate figures
    print("\n--- Generating Figures ---\n")

    if mode_str=="HISTORY MATCHING":
        # These plots only make sense for history matching runs
        print("1. Phi comparison plots...")
        plot_phi_comparison(data, output_dir)
        plot_phi_boxplot(data, output_dir)

        print("\n2. Time series observation plots...")
        plot_timeseries_obs(data, output_dir)

        # print("\n3. Parameter histogram plots...")
        # plot_parameter_histograms(data, output_dir)

        print("\n4. Parameter array statistics (posterior)...")
        plot_pname_array_statistics(data, 'pt_pe_fn', output_dir, 'ghb-conf-cond-fngr', 'ghb-conf-cond-fngr')
        plot_pname_array_statistics(data, 'pt_pe_fn', output_dir, 'ghb-conf-head-fngr', 'ghb-conf-head-fngr')

        print("\n5. Observation array statistics (posterior)...")
        # plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'confq', 'prior_confq')
        plot_oname_array_statistics(data, 'pt_oe_fn', output_dir, 'confq', 'posterior_confq')
        plot_oname_array_statistics(data, 'pt_oe_fn', output_dir, 'awq', 'posterior_awq')
        plot_oname_array_statistics(data, 'pt_oe_fn', output_dir, 'h-lyr', 'posterior_heads')

        print("\n6. 1:1 observed vs simulated plot...")
        plot_1to1(data, output_dir)

        print("\n7. 1:1 boxplot and scatter with uncertainty...")
        plot_1to1_boxplots(data, output_dir)
    else:
        # Prior-only mode: only generate array statistics
        print("1. Parameter array statistics (prior)...")
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'ghb-conf-cond-gr', 'ghb-conf-cond-gr')
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'ghb-conf-head-gr', 'ghb-conf-head-gr')
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'npfklayer1-fngr', 'npfklayer1-fngr', array_param=True)
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'npfklayer1-pp', 'npfklayer1-pp', array_param=True)
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'rch-fngr', 'rch-fngr', array_param=True)
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'rch-pp', 'rch-pp', array_param=True)
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'stosslayer1-pp', 'stosslayer1-pp', array_param=True)
        # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'stosslayer1-fngr', 'stosslayer1-fngr', array_param=True)

        # print("\n2. Observation array statistics (prior)...")
        # plot_heads_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'h-lyr', 'prior_heads')
        # plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'confq', 'prior_confq')
        # plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'awq', 'prior_awq')

        print("\n3. Check for prior data conflicts...")
        plot_pdc_budget_confined(data, output_dir)
        

    print(f"\n{'='*70}")
    print(f"POST-PROCESSING COMPLETE!")
    print(f"Figures saved to: {output_dir}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
