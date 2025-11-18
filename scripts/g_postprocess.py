"""
Post-processing script to generate figures from PEST++ IES results.
Adapted from notebooks/figure.ipynb

Usage:
    # Use default model from setup.py, auto-detect last iteration, all realizations:
    python scripts/g_postprocess.py

    # Specify model and iteration:
    python scripts/g_postprocess.py --model local_run34 --last-iter 8

    # Subset by phi percentile (default [10, 90]):
    python scripts/g_postprocess.py --subset-method phi --phi-percentile 10 90

    # Subset by mean head threshold:
    python scripts/g_postprocess.py --subset-method head --head-threshold 300

    # Use all realizations (no subsetting):
    python scripts/g_postprocess.py --subset-method all
"""

import os
import sys
import argparse
import pyemu
from g_helpers_run4 import *
from setup import MODEL_NAME


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate post-processing figures from PEST++ IES results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    parser.add_argument(
        '--model', '-m',
        type=str,
        default=MODEL_NAME,
        help=f'Model/run name (default: {MODEL_NAME} from setup.py)'
    )

    parser.add_argument(
        '--last-iter', '-i',
        type=int,
        default=None,
        help='Last iteration to use (default: auto-detect from PST file)'
    )

    parser.add_argument(
        '--subset-method',
        type=str,
        choices=['all', 'phi', 'head'],
        default='all',
        help='Method for subsetting realizations: all (no subset), phi (percentile), head (threshold)'
    )

    parser.add_argument(
        '--phi-percentile',
        type=float,
        nargs=2,
        default=[10, 90],
        metavar=('LOW', 'HIGH'),
        help='Phi percentile range for subsetting (default: 10 90)'
    )

    parser.add_argument(
        '--head-threshold',
        type=float,
        default=None,
        help='Mean head threshold for subsetting realizations'
    )

    return parser.parse_args()


def main():
    """Main post-processing function."""

    # Parse command-line arguments
    args = parse_args()
    run_name = args.model

    print(f"\n{'='*70}")
    print(f"POST-PROCESSING: {run_name}")
    print(f"{'='*70}\n")

    # Setup paths
    paths = setup_paths(run_name, m_d_str='master_ies', m_d_prior='master_ies')

    check_dir = paths['m_d']
    if not os.path.exists(check_dir):
        print(f"ERROR: Master IES directory not found: {check_dir}")
        print("Please run PEST++ IES first (e_pest_IES_HM.py)")
        return

    # Auto-detect last iteration if not specified
    if args.last_iter is None:
        pst_file = os.path.join(check_dir, f"{run_name}.pst")
        if os.path.exists(pst_file):
            pst = pyemu.Pst(pst_file)
            last_iter = pst.control_data.noptmax
            print(f"Auto-detected last iteration: {last_iter}")
        else:
            print(f"ERROR: PST file not found: {pst_file}")
            return
    else:
        last_iter = args.last_iter
        print(f"Using specified last iteration: {last_iter}")

    # Load data
    print("\nLoading data...")
    data = load_data(paths, run_name, last_iter=last_iter, use_prior_only=False, restart=False)

    output_dir = paths['output']

    # Handle realization subsetting based on method
    subset_oe = None
    subset_pe = None

    if args.subset_method == 'phi':
        print(f"\nSubsetting by phi percentile: {args.phi_percentile}")
        subset_indexes, phis = subset_realizations_by_phi(
            data,
            percentile=args.phi_percentile,
            from_iter=last_iter
        )
        suffix = f"_phi{int(args.phi_percentile[0])}-{int(args.phi_percentile[1])}"

        subset_oe_fn = os.path.join(output_dir, f'subset_oe{suffix}.csv')
        subset_pe_fn = os.path.join(output_dir, f'subset_pe{suffix}.csv')

        if not os.path.exists(subset_oe_fn):
            subset_oe, subset_pe = get_subset_ensembles(
                paths, data, run_name, subset_indexes,
                head_filter=False, output_suffix=suffix
            )
        else:
            print(f"Loading existing subset from {subset_oe_fn}")
            subset_oe = pd.read_csv(subset_oe_fn)
            subset_pe = pd.read_csv(subset_pe_fn)

    elif args.subset_method == 'head':
        if args.head_threshold is None:
            print("ERROR: --head-threshold required when using --subset-method head")
            return

        print(f"\nSubsetting by mean head threshold: {args.head_threshold}")
        subset_indexes, means = subset_realizations_by_mean_heads(
            data,
            'arr-h',
            min_threshold=args.head_threshold,
            from_iter=last_iter
        )
        suffix = f"_head{int(args.head_threshold)}"

        subset_oe_fn = os.path.join(output_dir, f'subset_oe{suffix}.csv')
        subset_pe_fn = os.path.join(output_dir, f'subset_pe{suffix}.csv')

        if not os.path.exists(subset_oe_fn):
            subset_oe, subset_pe = get_subset_ensembles(
                paths, data, run_name, subset_indexes,
                head_filter=True, output_suffix=suffix
            )
        else:
            print(f"Loading existing subset from {subset_oe_fn}")
            subset_oe = pd.read_csv(subset_oe_fn)
            subset_pe = pd.read_csv(subset_pe_fn)

    else:  # 'all'
        print("\nUsing all realizations (no subsetting)")

    # Generate plots with subset if specified
    if subset_oe is not None:
        plot_oname_array_statistics(data, subset_oe, output_dir, 'arr-h')

    # Generate figures
    print("\n--- Generating Figures ---\n")

    print_statistics(
        data,
        subset_oe=subset_oe,
        output=output_dir,
        oname='arr-h')

    # These plots only make sense for history matching runs
    print("1. Phi comparison plots...")
    plot_phi_comparison(data, output_dir)
    plot_phi_boxplot(data, output_dir)
    plot_phi_pie(data, output_dir)

    # plot parameter distributions
    print("\n2. Parameter distributions...")
    plot_parameter_distributions_subset(
        data,
        subset_pe=subset_pe,
        output_dir=output_dir)

    print("\n3. Prediction plots...")
    plot_histograms_spring_flux(
        data,
        subset_oe=subset_oe,
        output=output_dir,
        col_startswith='ts-sprflux')
    plot_timeseries_spring_flux(
        data,
        subset_oe=subset_oe,
        output=output_dir,
        col_startswith='ts-sprflux')

    print("\n4. Observation plots...")
    plot_timeseries_obs(
        data,
        subset_oe=subset_oe,
        output=output_dir,
        truth_file='truth/output.sample_heads.truth.csv',
        col_startswith='ts-heads')

    plot_budgets(
        data,
        subset_oe=subset_oe,
        output=output_dir,
        col_startswith='budget')

    print("\n5. Parameter histogram plots...")
    plot_parameter_histograms(
        data,
        subset_pe=subset_pe,
        output=output_dir)

    print("\n6. Observation array statistics (posterior)...")
    plot_oname_array_statistics(
        data,
        subset_oe=subset_oe,
        output=output_dir,
        oname_prefix='arr-h')
    plot_oname_array_statistics(data, output=output_dir, oname_prefix='arr-confq')
    plot_oname_array_statistics(data, output=output_dir, oname_prefix='awq')

    print("\n7. 1:1 observed vs simulated plot...")
    plot_1to1(data, output_dir)

    print("\n8. 1:1 boxplot and scatter with uncertainty...")
    plot_1to1_boxplots(data, output_dir)





    # Optional: Prior analysis plots (uncomment if needed)
    # Note: Use g_postprocess_prior.py for comprehensive prior analysis
    # print("\n9. Parameter array statistics (prior)...")
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'ghb-conf-cond-gr', 'ghb-conf-cond-gr')
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'ghb-conf-head-gr', 'ghb-conf-head-gr')
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'npfklayer1-fngr', 'npfklayer1-fngr', array_param=True)
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'npfklayer1-pp', 'npfklayer1-pp', array_param=True)
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'rch-fngr', 'rch-fngr', array_param=True)
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'rch-pp', 'rch-pp', array_param=True)
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'stosslayer1-pp', 'stosslayer1-pp', array_param=True)
    # plot_pname_array_statistics(data, 'pr_pe_fn', output_dir, 'stosslayer1-fngr', 'stosslayer1-fngr', array_param=True)
    #
    # print("\n10. Observation array statistics (prior)...")
    # plot_heads_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'h-lyr', 'prior_heads')
    # plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'confq', 'prior_confq')
    # plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'awq', 'prior_awq')
    #
    # print("\n11. Check for prior data conflicts...")
    # plot_pdc_budget_confined(data, output_dir)

    print(f"\n{'='*70}")
    print(f"POST-PROCESSING COMPLETE!")
    print(f"{'='*70}")
    print(f"Model: {run_name}")
    print(f"Last iteration: {last_iter}")
    if subset_oe is not None:
        print(f"Subset method: {args.subset_method}")
        if args.subset_method == 'phi':
            print(f"Phi percentile: {args.phi_percentile}")
        elif args.subset_method == 'head':
            print(f"Head threshold: {args.head_threshold}")
        print(f"Subset size: {len(subset_oe.columns)} realizations")
    else:
        print(f"Using all realizations (no subsetting)")
    print(f"Figures saved to: {output_dir}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
