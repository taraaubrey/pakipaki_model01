"""
pp_budget_boxplots.py - Create boxplots comparing prior vs posterior budget distributions

Creates boxplot comparisons for each budget type (Awanui, Poukawa, Spring, Confined,
Recharge, Storage) showing prior and posterior distributions for each stress period.

Usage:
    python scripts/pp_budget_boxplots.py run_name [options]

Examples:
    python scripts/pp_budget_boxplots.py local_run34 --post-iter 19
    python scripts/pp_budget_boxplots.py local_run34 --post-iter 19 --filter-file models/local_run34/filtered_realizations.csv
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyemu

# Add scripts directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import default model name from setup
try:
    from setup import MODEL_NAME as DEFAULT_MODEL_NAME
except ImportError:
    DEFAULT_MODEL_NAME = None


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Create boxplots comparing prior vs posterior budget distributions'
    )
    parser.add_argument(
        'model_name',
        type=str,
        nargs='?',
        default=DEFAULT_MODEL_NAME,
        help=f'Model name for file prefixes (default: {DEFAULT_MODEL_NAME} from setup.py)'
    )
    parser.add_argument('--run-name', '-r', type=str, default=None,
                       help='Run name for directory path (default: same as model_name)')
    parser.add_argument('--post-iter', type=int, default=19,
                       help='Posterior iteration (default: 19)')
    parser.add_argument('--filter-file', type=str, default=None,
                       help='Path to file with filtered realization indices')
    parser.add_argument('--suffix', type=str, default=None,
                       help='Suffix for output files (e.g., "_filtered")')
    return parser.parse_args()


def load_filtered_realizations(filter_file):
    """Load filtered realization indices from file."""
    if filter_file is None or not os.path.exists(filter_file):
        return None
    df = pd.read_csv(filter_file, comment='#')
    return df['realization'].tolist()


def get_budget_obs(obs_data):
    """
    Get budget observation names organized by type and kper.

    Supports both formats:
    - Old: oname='budget' with usecol patterns
    - New: oname='budget-sw', 'budget-rch', etc.

    Returns dict: {budget_type: {kper: obs_name}}
    """
    result = {}

    # Define budget types and their patterns
    budget_mappings = {
        'Awanui': {'old_patterns': ['awanui'], 'new_oname': 'budget-aw'},
        'Poukawa': {'old_patterns': ['poukawa'], 'new_oname': 'budget-pw'},
        'Spring': {'old_patterns': ['spring'], 'new_oname': 'budget-spring'},
        'Confined': {'old_patterns': ['confined'], 'new_oname': 'budget-conf'},
        'Recharge': {'old_patterns': ['recharge', 'rch'], 'new_oname': 'budget-rch'},
        'Surface Water': {'old_patterns': ['sw'], 'new_oname': 'budget-sw'},
        'Storage': {'old_patterns': ['stoss', 'storage'], 'new_oname': 'budget-storage'},
    }

    for btype, mapping in budget_mappings.items():
        result[btype] = {}

        # Try new format first (oname = 'budget-sw', 'budget-rch', etc.)
        new_oname = mapping['new_oname']
        new_obs = obs_data[obs_data['oname'] == new_oname].copy()

        if len(new_obs) > 0:
            new_obs['kper'] = pd.to_numeric(new_obs['kper'], errors='coerce')
            for idx, row in new_obs.iterrows():
                kper = int(row['kper']) if pd.notna(row['kper']) else 1
                result[btype][kper] = idx
        else:
            # Fall back to old format (oname = 'budget' with usecol)
            budget_obs = obs_data[obs_data['oname'] == 'budget'].copy()
            if len(budget_obs) > 0:
                budget_obs['kper'] = pd.to_numeric(budget_obs['kper'], errors='coerce')

                for idx, row in budget_obs.iterrows():
                    usecol = str(row.get('usecol', '')).lower()

                    # Check if this observation matches the budget type
                    for pattern in mapping['old_patterns']:
                        if pattern == usecol:
                            kper = int(row['kper']) if pd.notna(row['kper']) else 1
                            result[btype][kper] = idx
                            break

    # Remove empty budget types
    result = {k: v for k, v in result.items() if len(v) > 0}

    if len(result) == 0:
        print("  No budget observations found")

    return result


def load_truth_budgets():
    """
    Load truth budget values.

    Returns dict: {budget_type: {kper: value}}
    """
    truth_file = os.path.join('truth', 'output.budget.truth.csv')

    if not os.path.exists(truth_file):
        print(f"  Truth file not found: {truth_file}")
        return {}

    try:
        truth_df = pd.read_csv(truth_file)

        # Map column names to our budget type names
        col_mapping = {
            'awanui': 'Awanui',
            'poukawa': 'Poukawa',
            'spring': 'Spring',
            'confined': 'Confined',
            'rch': 'Recharge',
            'stoss': 'Storage'
        }

        # Convert to dict {budget_type: {kper: value}}
        truth_dict = {}
        for col, btype in col_mapping.items():
            if col in truth_df.columns:
                truth_dict[btype] = {}
                for _, row in truth_df.iterrows():
                    kper = int(row['kper'])
                    val = row[col]
                    truth_dict[btype][kper] = float(val)

        return truth_dict
    except Exception as e:
        print(f"  Error loading truth file: {e}")
        return {}


def create_budget_boxplots(run_name, model_name, post_iter=19, filter_file=None, suffix=None):
    """
    Create budget boxplot comparisons for prior vs posterior.
    """
    file_suffix = suffix if suffix else ''
    print(f"\n{'='*80}")
    print(f"BUDGET BOXPLOT ANALYSIS: {run_name}")
    print(f"{'='*80}")

    # Setup paths
    model_ws = os.path.join('models', run_name)
    master_dir = os.path.join(model_ws, 'pest', 'master_ies')
    output_dir = os.path.join(model_ws, 'pp_budget_analysis')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Load PST
    pst_file = os.path.join(master_dir, f'{model_name}.pst')
    pst = pyemu.Pst(pst_file)
    obs_data = pst.observation_data

    # Get budget observations
    budget_obs = get_budget_obs(obs_data)

    print("\nBudget observations found:")
    all_obs_names = []
    for btype, kper_dict in budget_obs.items():
        if kper_dict:
            print(f"  {btype}: {kper_dict}")
            all_obs_names.extend(kper_dict.values())

    if len(all_obs_names) == 0:
        print("Error: No budget observations found")
        return

    # Load truth values
    print("\nLoading truth values...")
    truth_values = load_truth_budgets()

    # Load filtered realizations
    filtered_ids = load_filtered_realizations(filter_file)
    if filtered_ids:
        print(f"\nUsing {len(filtered_ids)} filtered realizations")

    # Load observation ensembles
    print(f"\nLoading observation ensembles...")
    prior_file = os.path.join(master_dir, f'{model_name}.0.obs.csv')
    post_file = os.path.join(master_dir, f'{model_name}.{post_iter}.obs.csv')

    # Read only needed columns
    print(f"  Reading prior (iter 0)...")
    prior_header = pd.read_csv(prior_file, nrows=0)
    cols_to_read = [prior_header.columns[0]]
    cols_to_read.extend([col for col in all_obs_names if col in prior_header.columns])
    prior_oe = pd.read_csv(prior_file, usecols=cols_to_read, index_col=0)

    print(f"  Reading posterior (iter {post_iter})...")
    post_oe = pd.read_csv(post_file, usecols=cols_to_read, index_col=0)

    # Filter realizations
    if filtered_ids:
        # Convert index to same type as filtered_ids for consistent comparison
        # Handle mixed types by converting to strings then back
        prior_index = [str(i) for i in prior_oe.index]
        post_index = [str(i) for i in post_oe.index]
        filtered_str = [str(r) for r in filtered_ids]

        prior_mask = [i for i, idx in enumerate(prior_index) if idx in filtered_str]
        post_mask = [i for i, idx in enumerate(post_index) if idx in filtered_str]

        prior_oe = prior_oe.iloc[prior_mask]
        post_oe = post_oe.iloc[post_mask]
        print(f"  Filtered to {len(prior_oe)} prior, {len(post_oe)} posterior realizations")

    # Get budget types that have data
    budget_types_with_data = [bt for bt, kd in budget_obs.items() if len(kd) > 0]

    # Create boxplots for each budget type
    print("\nCreating boxplots...")

    for btype in budget_types_with_data:
        kper_dict = budget_obs[btype]
        kpers = sorted(kper_dict.keys())

        if len(kpers) == 0:
            continue

        # Create figure
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))

        # Prepare data for boxplot
        box_data = []
        positions = []
        colors = []
        labels = []

        pos = 0
        for kper in kpers:
            obs_name = kper_dict[kper]

            if obs_name not in prior_oe.columns:
                continue

            # Prior data
            prior_vals = prior_oe[obs_name].values
            box_data.append(prior_vals)
            positions.append(pos)
            colors.append('#1f77b4')
            labels.append(f'kper {kper}\nPrior')

            # Posterior data
            post_vals = post_oe[obs_name].values
            box_data.append(post_vals)
            positions.append(pos + 0.4)
            colors.append('#ff7f0e')
            labels.append(f'kper {kper}\nPost')

            pos += 1

        # Create boxplot
        bp = ax.boxplot(box_data, positions=positions, widths=0.35, patch_artist=True,
                       showfliers=False)

        # Color the boxes
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)

        # Add truth values
        if btype in truth_values:
            for i, kper in enumerate(kpers):
                if kper in truth_values[btype]:
                    truth_val = truth_values[btype][kper]
                    # Draw horizontal line spanning both prior and posterior boxes
                    x_start = i - 0.25
                    x_end = i + 0.65
                    ax.hlines(truth_val, x_start, x_end, colors='red', linestyles='--',
                             linewidth=2, label='Truth' if i == 0 else '')

        # Customize plot
        ax.set_xlabel('Stress Period', fontsize=10)
        ax.set_ylabel('Budget (m³/d)', fontsize=10)
        ax.set_title(f'{btype} Budget - Prior vs Posterior', fontsize=12, fontweight='bold')

        # Set x-ticks at the center of each kper group
        tick_positions = [i + 0.2 for i in range(len(kpers))]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels([f'kper {k}' for k in kpers])

        ax.grid(True, alpha=0.3, axis='y')

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#1f77b4', alpha=0.7, label='Prior'),
            Patch(facecolor='#ff7f0e', alpha=0.7, label='Posterior'),
        ]
        if btype in truth_values:
            from matplotlib.lines import Line2D
            legend_elements.append(Line2D([0], [0], color='red', linestyle='--', linewidth=2, label='Truth'))

        ax.legend(handles=legend_elements, loc='best', fontsize=9)

        plt.tight_layout()

        # Save figure
        output_path = os.path.join(output_dir, f'{run_name}_{btype.lower()}_boxplot{file_suffix}.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {output_path}")

    # Create separate figures for prior and posterior
    print("\nCreating separate prior and posterior plots...")

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D

    # Define kper order and colors
    kper_order = [1, 3, 2, 4]
    kper_colors = {1: '#1f77b4', 2: '#1f77b4', 3: '#9467bd', 4: '#9467bd'}  # Purple for kper 3, 4

    # Panel labels
    panel_labels = ['A', 'B', 'C', 'D', 'E', 'F']

    # Collect summary statistics
    summary_rows = []

    for iteration_name, oe_data in [('prior', prior_oe), ('posterior', post_oe)]:
        n_types = len(budget_types_with_data)
        fig, axes = plt.subplots(2, 3, figsize=(7.28, 7.28/2))
        axes_flat = axes.flatten()

        for idx, btype in enumerate(budget_types_with_data):
            if idx >= 6:
                break

            ax = axes_flat[idx]
            kper_dict = budget_obs[btype]

            # Use custom kper order, filtering to only available kpers
            kpers = [k for k in kper_order if k in kper_dict]

            if len(kpers) == 0:
                ax.set_visible(False)
                continue

            # Collect data for boxplot and calculate IQR limits
            box_data = []
            all_vals = []

            for i, kper in enumerate(kpers):
                obs_name = kper_dict[kper]

                if obs_name not in oe_data.columns:
                    box_data.append([])
                    continue

                vals = oe_data[obs_name].values
                box_data.append(vals)
                all_vals.extend(vals)

                # Collect summary statistics
                summary_rows.append({
                    'iteration': iteration_name,
                    'budget_type': btype,
                    'kper': kper,
                    'n': len(vals),
                    'min': np.min(vals),
                    'p05': np.percentile(vals, 5),
                    'p25': np.percentile(vals, 25),
                    'median': np.median(vals),
                    'mean': np.mean(vals),
                    'p75': np.percentile(vals, 75),
                    'p95': np.percentile(vals, 95),
                    'max': np.max(vals),
                    'std': np.std(vals),
                })

            # Calculate IQR-based y-axis limits
            if len(all_vals) > 0:
                all_vals = np.array(all_vals)
                q25 = np.percentile(all_vals, 25)
                q75 = np.percentile(all_vals, 75)
                iqr = q75 - q25
                y_min = max(np.min(all_vals), q25 - 1.5 * iqr)
                y_max = min(np.max(all_vals), q75 + 1.5 * iqr)
                # Add small margin
                margin = (y_max - y_min) * 0.05
                ax.set_ylim(y_min - margin, y_max + margin)

            # Plot scatter dots first (very small, black)
            for i, kper in enumerate(kpers):
                obs_name = kper_dict[kper]

                if obs_name not in oe_data.columns:
                    continue

                vals = oe_data[obs_name].values
                # Add jitter to x positions for visibility
                x_jitter = np.random.normal(i, 0.08, len(vals))
                ax.scatter(x_jitter, vals, c='black', s=0.5, alpha=0.3, zorder=1)

            # Overlay boxplot
            bp = ax.boxplot(box_data, positions=range(len(kpers)), widths=0.5,
                           patch_artist=True, showfliers=False, zorder=2)

            # Color the boxes
            for i, (patch, kper) in enumerate(zip(bp['boxes'], kpers)):
                color = kper_colors.get(kper, '#1f77b4')
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            # Style whiskers and medians
            for element in ['whiskers', 'caps']:
                for line in bp[element]:
                    line.set_color('black')
                    line.set_linewidth(0.8)
            for median in bp['medians']:
                median.set_color('black')
                median.set_linewidth(1.2)

            # Customize plot
            ax.set_ylabel('Budget (m³/d)', fontsize=6)
            # Left-justified title with panel label
            ax.set_title(f'{panel_labels[idx]}: {btype}', fontsize=6, fontweight='bold', loc='left')

            ax.set_xticks(range(len(kpers)))
            ax.set_xticklabels([f'{k}' for k in kpers], fontsize=6)
            ax.set_xlabel('kper', fontsize=6)

            ax.grid(True, alpha=0.3, axis='y')
            ax.tick_params(labelsize=6)

        # Hide unused subplots
        for idx in range(len(budget_types_with_data), 6):
            axes_flat[idx].set_visible(False)

        # Add legend to figure
        legend_elements = [
            Patch(facecolor='#1f77b4', alpha=0.7, label='kper 1, 2'),
            Patch(facecolor='#9467bd', alpha=0.7, label='kper 3, 4'),
        ]
        fig.legend(handles=legend_elements, loc='upper right', fontsize=6, bbox_to_anchor=(0.98, 0.98))

        plt.tight_layout()

        output_path = os.path.join(output_dir, f'{run_name}_budget_boxplots_{iteration_name}{file_suffix}.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        print(f"  Saved: {output_path}")

    # Save summary statistics table
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_path = os.path.join(output_dir, f'{run_name}_budget_summary{file_suffix}.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"  Saved summary: {summary_path}")

    print(f"\nAll plots saved to: {output_dir}")


if __name__ == '__main__':
    args = parse_args()

    model_name = args.model_name
    run_name = args.run_name if args.run_name else model_name

    print(f"Model name: {model_name}")
    if run_name != model_name:
        print(f"Run name (directory): {run_name}")

    create_budget_boxplots(
        run_name=run_name,
        model_name=model_name,
        post_iter=args.post_iter,
        filter_file=args.filter_file,
        suffix=args.suffix
    )
