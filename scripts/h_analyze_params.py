"""
Analyze PEST++ IES parameter change summary (PCS) results.
Run this after e_pest_IES_HM.py to understand parameter evolution and behavior.

Usage:
    python scripts/h_analyze_params.py [model_name]

    If model_name is not provided, uses MODEL_NAME from setup.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import glob
import pyemu

# Add scripts directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

# Import setup to get model name
from setup import MODEL_NAME as DEFAULT_MODEL_NAME


def categorize_parameter(pname):
    """
    Categorize parameter by type based on pname prefix.

    Args:
        pname: Parameter name (e.g., 'npfklayer1-gr')

    Returns:
        tuple: (category, subcategory)
    """
    pname_lower = str(pname).lower()

    if 'npfk' in pname_lower:
        if '-gr' in pname_lower:
            return ('Hydraulic Conductivity', 'K-grid')
        elif '-pp' in pname_lower:
            return ('Hydraulic Conductivity', 'K-pilot')
        elif '-cn' in pname_lower:
            return ('Hydraulic Conductivity', 'K-constant')
        else:
            return ('Hydraulic Conductivity', 'K-other')

    elif 'stoss' in pname_lower or 'sto' in pname_lower:
        if '-gr' in pname_lower:
            return ('Storage', 'S-grid')
        elif '-pp' in pname_lower:
            return ('Storage', 'S-pilot')
        elif '-cn' in pname_lower:
            return ('Storage', 'S-constant')
        else:
            return ('Storage', 'S-other')

    elif 'rch' in pname_lower:
        if 'rchss' in pname_lower or 'ss' in pname_lower:
            return ('Recharge', 'RCH-steady')
        elif 'rchtr' in pname_lower or 'tr' in pname_lower:
            return ('Recharge', 'RCH-transient')
        else:
            return ('Recharge', 'RCH-other')

    elif 'ghb' in pname_lower:
        # Determine GHB type
        ghb_type = 'GHB-other'
        if 'ghbaw' in pname_lower:
            ghb_type = 'GHB-Awanui'
        elif 'ghbpw' in pname_lower:
            ghb_type = 'GHB-Poukawa'
        elif 'ghbspring' in pname_lower:
            ghb_type = 'GHB-Spring'
        elif 'ghbconf' in pname_lower:
            ghb_type = 'GHB-Confined'

        # Determine property
        if '-cond' in pname_lower:
            return ('GHB Conductance', ghb_type + '-cond')
        elif '-head' in pname_lower:
            return ('GHB Head', ghb_type + '-head')
        else:
            return ('GHB Other', ghb_type)

    else:
        return ('Other', 'Other')


def analyze_parameters(model_name):
    """
    Analyze parameter change summary results for a given model run.

    Args:
        model_name: Name of the model run (e.g., 'local_run33')
    """
    print("="*100)
    print(f"PEST++ IES PARAMETER ANALYSIS: {model_name}")
    print("="*100)

    # Set paths
    master_dir = os.path.join('models', model_name, 'pest', 'master_ies')
    output_dir = os.path.join('models', model_name)

    if not os.path.exists(master_dir):
        print(f"ERROR: Directory not found: {master_dir}")
        print("Make sure you've run e_pest_IES_HM.py first!")
        return

    # Find all PCS files
    pcs_pattern = os.path.join(master_dir, f'{model_name}.*.pcs.csv')
    pcs_files = sorted(glob.glob(pcs_pattern))

    if len(pcs_files) == 0:
        print(f"ERROR: No PCS files found matching: {pcs_pattern}")
        return

    print(f"\nFound {len(pcs_files)} PCS files")

    # Load parameter group mapping
    mapping_file = os.path.join(master_dir, 'par_group_mapping.csv')
    if os.path.exists(mapping_file):
        mapping = pd.read_csv(mapping_file)
        # Create dict from pargp to pname
        pargp_to_pname = dict(zip(mapping['pargp'], mapping['pname']))
        print(f"Loaded parameter group mapping: {len(pargp_to_pname)} groups")
    else:
        print(f"Mapping file not found, extracting from PST file...")
        # Load PST file to get parameter information
        pst_file = os.path.join(master_dir, f'{model_name}.pst')
        if os.path.exists(pst_file):
            pst = pyemu.Pst(pst_file)
            # Create mapping from pargp to pname using the pname column
            pargp_to_pname = {}
            for pargp in pst.par_groups:
                # Get all parameters in this group
                pars = pst.parameter_data[pst.parameter_data['pargp'] == pargp]
                if len(pars) > 0:
                    # Use the pname column directly (should be consistent within a pargp)
                    if 'pname' in pars.columns:
                        pname = pars['pname'].iloc[0]
                        pargp_to_pname[pargp] = pname
                    else:
                        # Fallback to parsing from index
                        first_pname = pars.index[0]
                        if 'pname:' in first_pname:
                            pname = first_pname.split('pname:')[1].split('_')[0]
                            pargp_to_pname[pargp] = pname
                        else:
                            pargp_to_pname[pargp] = pargp
            print(f"Extracted parameter group mapping from PST: {len(pargp_to_pname)} groups")
        else:
            print(f"ERROR: PST file not found: {pst_file}")
            pargp_to_pname = {}

    # Load all PCS files
    all_pcs = []
    for pcs_file in pcs_files:
        # Extract iteration number from filename
        # Format: {model_name}.{iteration}.pcs.csv or {model_name}.{iteration}.reinflate.pcs.csv
        basename = os.path.basename(pcs_file)
        # Remove .pcs.csv suffix, then remove .reinflate if present
        name_parts = basename.replace('.pcs.csv', '').replace('.reinflate', '').split('.')
        # Get the last numeric part
        iter_num = int(name_parts[-1])

        pcs = pd.read_csv(pcs_file)
        pcs['iteration'] = iter_num

        # Add pname from mapping
        pcs['pname'] = pcs['group'].map(pargp_to_pname)
        # If not in mapping, use group name
        pcs['pname'] = pcs['pname'].fillna(pcs['group'])

        # Categorize parameters
        pcs[['category', 'subcategory']] = pcs['pname'].apply(
            lambda x: pd.Series(categorize_parameter(x))
        )

        all_pcs.append(pcs)

    # Combine all iterations
    pcs_df = pd.concat(all_pcs, ignore_index=True)

    print(f"\nParameter categories found:")
    for cat in pcs_df['category'].unique():
        count = pcs_df[pcs_df['category'] == cat]['subcategory'].nunique()
        print(f"  {cat}: {count} subcategories")

    # Load phi data for comparison
    phi_actual_file = os.path.join(master_dir, f'{model_name}.phi.actual.csv')
    if os.path.exists(phi_actual_file):
        phi_actual = pd.read_csv(phi_actual_file)
        print(f"\nLoaded phi data: {len(phi_actual)} iterations")
    else:
        print(f"\nWARNING: Phi file not found: {phi_actual_file}")
        phi_actual = None

    # Calculate statistics by parameter group across iterations
    print("\n" + "="*100)
    print("PARAMETER CHANGE SUMMARY BY ITERATION")
    print("="*100)

    summary_stats = []
    for _, row in pcs_df.iterrows():
        # Get initial and current values
        initial_cv = row['initial_cv']
        current_cv = row['current_cv']
        cv_reduction = ((initial_cv - current_cv) / initial_cv * 100) if initial_cv > 0 else 0

        summary_stats.append({
            'Iteration': row['iteration'],
            'Category': row['category'],
            'Subcategory': row['subcategory'],
            'Group': row['group'],
            'Pname': row['pname'],
            'Mean_Change': row['mean_change'],
            'Std_Change': row['std_change'],
            'Initial_CV': initial_cv,
            'Current_CV': current_cv,
            'CV_Reduction_%': cv_reduction,
            'Pct_at_Lower_Bound': row['percent_at_near_lbound'],
            'Pct_at_Upper_Bound': row['percent_at_near_ubound'],
        })

    stats_df = pd.DataFrame(summary_stats)

    # Get final iteration stats
    final_iter = stats_df['Iteration'].max()
    final_stats = stats_df[stats_df['Iteration'] == final_iter].copy()
    final_stats = final_stats.sort_values('CV_Reduction_%', ascending=False)

    print("\nFinal Iteration Parameter Statistics (sorted by CV reduction):")
    print(final_stats[['Category', 'Subcategory', 'Pname', 'Mean_Change',
                       'CV_Reduction_%', 'Pct_at_Lower_Bound', 'Pct_at_Upper_Bound']].to_string(index=False))

    # Analyze by category
    print("\n" + "="*100)
    print("PARAMETER EVOLUTION BY CATEGORY")
    print("="*100)

    for category in stats_df['Category'].unique():
        cat_df = stats_df[stats_df['Category'] == category]
        print(f"\n{category}:")

        for subcategory in cat_df['Subcategory'].unique():
            subcat_df = cat_df[cat_df['Subcategory'] == subcategory]

            # Get final iteration row for this subcategory
            final_row = subcat_df[subcat_df['Iteration'] == final_iter].iloc[0] if len(subcat_df[subcat_df['Iteration'] == final_iter]) > 0 else None

            if final_row is not None:
                initial_cv = final_row['Initial_CV']
                current_cv = final_row['Current_CV']
                cv_reduction_pct = final_row['CV_Reduction_%']

                print(f"  {subcategory}:")
                print(f"    Initial CV: {initial_cv:.2f}")
                print(f"    Final CV: {current_cv:.2f}")
                print(f"    CV Reduction: {cv_reduction_pct:.1f}%")
                print(f"    Final Mean Change: {final_row['Mean_Change']:.2f}")
                print(f"    At Lower Bound: {final_row['Pct_at_Lower_Bound']:.1f}%")
                print(f"    At Upper Bound: {final_row['Pct_at_Upper_Bound']:.1f}%")

    # Create visualization
    print("\n" + "="*100)
    print("CREATING VISUALIZATIONS...")
    print("="*100)

    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Color scheme for parameter categories
    category_colors = {
        'Hydraulic Conductivity': '#e74c3c',
        'Storage': '#3498db',
        'Recharge': '#2ecc71',
        'GHB Conductance': '#f39c12',
        'GHB Head': '#9b59b6',
        'GHB Other': '#95a5a6',
        'Other': '#34495e'
    }

    # Plot 1: Phi evolution with parameter change overlay
    if phi_actual is not None:
        ax1 = fig.add_subplot(gs[0, :])

        # Phi on left axis
        ax1.plot(phi_actual['iteration'], phi_actual['mean'], 'o-', linewidth=3,
                markersize=8, color='black', label='Mean Phi', zorder=3)
        ax1.fill_between(phi_actual['iteration'],
                        phi_actual['mean'] - phi_actual['standard_deviation'],
                        phi_actual['mean'] + phi_actual['standard_deviation'],
                        alpha=0.2, color='black')
        ax1.set_ylabel('Total Phi', fontsize=12, fontweight='bold')
        ax1.set_yscale('log')
        ax1.set_xlabel('Iteration', fontsize=12, fontweight='bold')
        ax1.legend(loc='upper left', fontsize=10)
        ax1.grid(True, alpha=0.3)

        # Average parameter change on right axis
        ax1b = ax1.twinx()
        avg_change = stats_df.groupby('Iteration')['Mean_Change'].mean()
        ax1b.plot(avg_change.index, avg_change.values, 's-', linewidth=2,
                 markersize=6, color='#e74c3c', label='Avg Parameter Change', alpha=0.7)
        ax1b.set_ylabel('Average Parameter Change', fontsize=12, fontweight='bold', color='#e74c3c')
        ax1b.tick_params(axis='y', labelcolor='#e74c3c')
        ax1b.legend(loc='upper right', fontsize=10)

        ax1.set_title(f'Phi Evolution vs Parameter Changes: {model_name}',
                     fontsize=14, fontweight='bold')

    # Plot 2: CV reduction by parameter category
    ax2 = fig.add_subplot(gs[1, 0])
    for category in stats_df['Category'].unique():
        cat_data = stats_df[stats_df['Category'] == category]
        # Average CV across subcategories
        cv_by_iter = cat_data.groupby('Iteration')['Current_CV'].mean()
        color = category_colors.get(category, '#95a5a6')
        ax2.plot(cv_by_iter.index, cv_by_iter.values, 'o-', linewidth=2,
                label=category, color=color, markersize=5)

    ax2.set_xlabel('Iteration', fontsize=11)
    ax2.set_ylabel('Coefficient of Variation', fontsize=11)
    ax2.set_title('Parameter Uncertainty (CV) Evolution', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=8, loc='best')
    ax2.grid(True, alpha=0.3)

    # Plot 3: Mean parameter change by category
    ax3 = fig.add_subplot(gs[1, 1])
    for category in stats_df['Category'].unique():
        cat_data = stats_df[stats_df['Category'] == category]
        change_by_iter = cat_data.groupby('Iteration')['Mean_Change'].mean()
        color = category_colors.get(category, '#95a5a6')
        ax3.plot(change_by_iter.index, change_by_iter.values, 'o-', linewidth=2,
                label=category, color=color, markersize=5)

    ax3.set_xlabel('Iteration', fontsize=11)
    ax3.set_ylabel('Mean Change Magnitude', fontsize=11)
    ax3.set_title('Parameter Change Magnitude', fontsize=12, fontweight='bold')
    ax3.legend(fontsize=8, loc='best')
    ax3.grid(True, alpha=0.3)

    # Plot 4: Boundary hitting analysis
    ax4 = fig.add_subplot(gs[1, 2])
    final_bounds = final_stats.groupby('Category').agg({
        'Pct_at_Lower_Bound': 'mean',
        'Pct_at_Upper_Bound': 'mean'
    })

    x = np.arange(len(final_bounds))
    width = 0.35
    categories_list = list(final_bounds.index)
    colors_list = [category_colors.get(cat, '#95a5a6') for cat in categories_list]

    bars1 = ax4.bar(x - width/2, final_bounds['Pct_at_Lower_Bound'], width,
                   label='At Lower Bound', color=colors_list, alpha=0.6, edgecolor='black')
    bars2 = ax4.bar(x + width/2, final_bounds['Pct_at_Upper_Bound'], width,
                   label='At Upper Bound', color=colors_list, alpha=0.9, edgecolor='black')

    ax4.set_ylabel('Percentage (%)', fontsize=11)
    ax4.set_title('Parameters at Bounds (Final Iteration)', fontsize=12, fontweight='bold')
    ax4.set_xticks(x)
    ax4.set_xticklabels(categories_list, rotation=45, ha='right', fontsize=9)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3, axis='y')
    ax4.axhline(y=5, color='orange', linestyle='--', linewidth=1, alpha=0.5, label='5% threshold')

    # Plot 5: CV reduction percentage by subcategory (final iteration)
    ax5 = fig.add_subplot(gs[2, :])
    final_sorted = final_stats.sort_values('CV_Reduction_%', ascending=False)

    bar_colors = [category_colors.get(cat, '#95a5a6') for cat in final_sorted['Category']]
    bars = ax5.bar(range(len(final_sorted)), final_sorted['CV_Reduction_%'],
                  color=bar_colors, alpha=0.8, edgecolor='black')

    ax5.set_xticks(range(len(final_sorted)))
    ax5.set_xticklabels(final_sorted['Subcategory'], rotation=45, ha='right', fontsize=9)
    ax5.set_ylabel('CV Reduction (%)', fontsize=11)
    ax5.set_title('Parameter Uncertainty Reduction from Initial (Final Iteration)',
                 fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3, axis='y')
    ax5.axhline(y=50, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    ax5.axhline(y=25, color='gray', linestyle=':', linewidth=1, alpha=0.5)

    # Add legend for categories
    from matplotlib.patches import Rectangle
    legend_elements = [Rectangle((0,0),1,1, fc=color, ec='black', label=cat)
                      for cat, color in category_colors.items()
                      if cat in final_sorted['Category'].values]
    ax5.legend(handles=legend_elements, loc='upper right', fontsize=9, ncol=2)

    # Save figure
    fig_path = os.path.join(output_dir, f'{model_name}_parameter_analysis.png')
    plt.savefig(fig_path, dpi=200, bbox_inches='tight')
    print(f"Figure saved to: {fig_path}")
    plt.close()

    # Save summary to CSV
    csv_path = os.path.join(output_dir, f'{model_name}_parameter_summary.csv')
    stats_df.to_csv(csv_path, index=False)
    print(f"Summary table saved to: {csv_path}")

    # Create markdown summary
    md_path = os.path.join(output_dir, f'{model_name}_PARAMETER_ANALYSIS.md')
    with open(md_path, 'w') as f:
        f.write(f"# Parameter Analysis Summary: {model_name}\n\n")
        f.write("## Overview\n\n")
        f.write(f"Total iterations completed: {final_iter}\n\n")
        f.write(f"Total parameter groups analyzed: {len(stats_df['Group'].unique())}\n\n")
        f.write(f"Parameter categories: {', '.join(stats_df['Category'].unique())}\n\n")

        if phi_actual is not None:
            f.write(f"Initial phi (mean): {phi_actual.loc[0, 'mean']:.2e}\n\n")
            f.write(f"Final phi (mean): {phi_actual.loc[phi_actual.index[-1], 'mean']:.2e}\n\n")
            phi_reduction = ((phi_actual.loc[0, 'mean'] - phi_actual.loc[phi_actual.index[-1], 'mean']) /
                           phi_actual.loc[0, 'mean'] * 100)
            f.write(f"Phi reduction: {phi_reduction:.2f}%\n\n")

        f.write("## Parameter Evolution by Category\n\n")
        for category in stats_df['Category'].unique():
            f.write(f"### {category}\n\n")
            cat_final = final_stats[final_stats['Category'] == category]
            if len(cat_final) > 0:
                f.write(cat_final[['Subcategory', 'Pname', 'Mean_Change', 'CV_Reduction_%',
                                  'Pct_at_Lower_Bound', 'Pct_at_Upper_Bound']].to_markdown(index=False))
                f.write("\n\n")

        f.write("## Boundary Hitting Diagnostics\n\n")
        f.write("Parameters at bounds may indicate over-parameterization or improper bounds.\n\n")

        high_lower = final_stats[final_stats['Pct_at_Lower_Bound'] > 5]
        high_upper = final_stats[final_stats['Pct_at_Upper_Bound'] > 5]

        if len(high_lower) > 0:
            f.write("### Groups with >5% at Lower Bound:\n\n")
            f.write(high_lower[['Subcategory', 'Pname', 'Pct_at_Lower_Bound']].to_markdown(index=False))
            f.write("\n\n")

        if len(high_upper) > 0:
            f.write("### Groups with >5% at Upper Bound:\n\n")
            f.write(high_upper[['Subcategory', 'Pname', 'Pct_at_Upper_Bound']].to_markdown(index=False))
            f.write("\n\n")

        f.write("## Uncertainty Reduction Summary\n\n")
        cv_summary = final_stats.groupby('Category')['CV_Reduction_%'].agg(['mean', 'min', 'max'])
        f.write(cv_summary.to_markdown())
        f.write("\n\n")

    print(f"Markdown summary saved to: {md_path}")

    print("\n" + "="*100)
    print("ANALYSIS COMPLETE!")
    print("="*100)

    return stats_df


if __name__ == '__main__':
    # Get model name from command line or use default
    if len(sys.argv) > 1:
        model_name = sys.argv[1]
    else:
        model_name = DEFAULT_MODEL_NAME
        print(f"No model name provided, using MODEL_NAME from setup.py: {model_name}")

    analyze_parameters(model_name)
