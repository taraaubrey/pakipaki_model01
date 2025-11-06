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

from g_helpers_prior import *
from setup import MODEL_NAME


def main():

    # Use MODEL_NAME from setup.py
    # run_name = MODEL_NAME
    run_name = 'local_run10'

    # args.prior_only = False

    print(f"\n{'='*70}")
    print(f"PRIOR POST-PROCESSING FIGURES FOR: {run_name}")
    print(f"{'='*70}\n")

    # Setup paths
    paths = setup_paths(run_name)

    # Load data
    print("Loading data...")
    data = load_prior_data(paths['m_d_pr'], run_name)
    output_dir = paths['output']


    # Generate figures
    print("\n--- Generating Figures ---\n")

    plot_oname_array_statistics(data, 'pr_oe_fn', output_dir, 'h-lyr', 'prior_heads')

    print("1. Prior data conflict...")
    plot_prior_data_conflict(data, output_dir)




    print(f"\n{'='*70}")
    print(f"POST-PROCESSING COMPLETE!")
    print(f"Figures saved to: {output_dir}")
    print(f"{'='*70}\n")

if __name__ == "__main__":
    main()
