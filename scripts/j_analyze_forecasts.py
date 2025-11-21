"""
PEST++ IES Forecast Analysis

Analyzes forecast/prediction uncertainty reduction and evolution across IES iterations.

Forecasts include water budget components for prediction periods:
- kper 1-4 (steady-state and transient periods)

Author: Claude Code
Date: 2025-11-19
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pyemu
from pathlib import Path

# Add scripts directory to path
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from setup import MODEL_NAME


def load_forecast_data(model_dir, model_name):
    """Load forecast data from PEST++ IES results.

    Args:
        model_dir: Path to model directory
        model_name: Name of the model/run

    Returns:
        dict: Dictionary containing forecast data
    """
    ies_dir = os.path.join(model_dir, 'pest', 'master_ies')
    prior_dir = os.path.join(model_dir, 'pest', 'master_ies_prior')

    # Load PST file
    pst_file = os.path.join(ies_dir, f'{model_name}.pst')
    if not os.path.exists(pst_file):
        raise FileNotFoundError(f"PST file not found: {pst_file}")

    pst = pyemu.Pst(pst_file)

    # Get forecast names
    forecasts = pst.forecast_names
    if len(forecasts) == 0:
        raise ValueError("No forecasts found in PST file")

    # Load observation data to get forecast metadata
    obs_data_file = os.path.join(ies_dir, f'{model_name}.obs_data.csv')
    obs_data = pd.read_csv(obs_data_file, index_col=0, low_memory=False)

    # Get forecast observation data
    forecast_obs = obs_data.loc[forecasts]

    # Find all observation ensemble files
    obs_files = sorted([f for f in os.listdir(ies_dir) if f.endswith('.obs.csv') and not 'noise' in f])
    iterations = []
    for f in obs_files:
        try:
            iter_num = int(f.replace(f'{model_name}.', '').replace('.obs.csv', ''))
            iterations.append(iter_num)
        except ValueError:
            continue

    iterations = sorted(iterations)

    # Load prior ensemble - only forecast columns for speed
    prior_file = os.path.join(prior_dir, f'{model_name}.0.obs.csv')
    if not os.path.exists(prior_file):
        prior_file = os.path.join(ies_dir, f'{model_name}.0.obs.csv')

    # Convert forecasts to list if needed and create column list for usecols
    forecast_cols = ['real_name'] + list(forecasts)
    prior_df = pd.read_csv(prior_file, index_col=0, usecols=forecast_cols)
    prior_forecasts = prior_df

    # Load posterior ensemble (last iteration) - only forecast columns
    last_iter = iterations[-1]
    post_file = os.path.join(ies_dir, f'{model_name}.{last_iter}.obs.csv')
    post_df = pd.read_csv(post_file, index_col=0, usecols=forecast_cols)
    post_forecasts = post_df

    # Load phi data to identify reinflation iterations
    phi_file = os.path.join(ies_dir, f'{model_name}.phi.actual.csv')
    phi_df = pd.read_csv(phi_file, index_col=0)

    # Detect reinflation iterations (where phi increases)
    mean_phi = phi_df.mean(axis=0)

    # Filter to only numeric iteration columns
    iter_cols = []
    for col in mean_phi.index:
        try:
            iter_cols.append(int(col))
        except (ValueError, TypeError):
            continue

    iter_cols = sorted(iter_cols)
    reinflate_iters = []

    for i in range(1, len(iter_cols)):
        curr_iter = iter_cols[i]
        prev_iter = iter_cols[i-1]
        if str(curr_iter) in mean_phi.index and str(prev_iter) in mean_phi.index:
            if mean_phi[str(curr_iter)] > mean_phi[str(prev_iter)]:
                reinflate_iters.append(curr_iter)

    return {
        'pst': pst,
        'forecasts': forecasts,
        'forecast_obs': forecast_obs,
        'iterations': iterations,
        'reinflate_iters': reinflate_iters,
        'prior': prior_forecasts,
        'posterior': post_forecasts,
        'ies_dir': ies_dir,
        'model_name': model_name,
        'last_iter': last_iter
    }


def categorize_forecast(row):
    """Categorize forecast by type.

    Args:
        row: Row from forecast_obs dataframe

    Returns:
        tuple: (category, subcategory, component)
    """
    longname = str(row['longname']).lower()
    obgnme = str(row['obgnme']).lower()
    oname = str(row['oname']).lower()
    kper = row['kper']
    kstp = row['kstp']
    i = row['i'] if 'i' in row.index else -1
    j = row['j'] if 'j' in row.index else -1

    # Determine subcategory based on stress period and timestep
    if kper == 1:
        subcategory = 'kper1'
        period = 'kper1'
    elif kper == 2:
        subcategory = 'kper2'
        period = 'kper2'
    elif kper == 3:
        if kstp <= 1:
            subcategory = 'kper3-ss'
            period = 'kper3'
        else:
            subcategory = f'kper3-ts{kstp}'
            period = f'kper3-TS{kstp}'
    elif kper == 4:
        if kstp <= 1:
            subcategory = 'kper4-ss'
            period = 'kper4'
        else:
            subcategory = f'kper4-ts{kstp}'
            period = f'kper4-TS{kstp}'
    else:
        subcategory = f'kper{kper}'
        period = f'kper{kper}'

    # Extract component name based on oname first, then longname
    if oname == 'ts-heads':
        # Timeseries heads forecasts
        component = f'Head TS (i{i:.0f},j{j:.0f})'
        category = 'Forecast-Heads'
    elif oname == 'ts-sprflux':
        # Spring flux timeseries (e.g., i=36, j=32)
        component = f'Spring Flux (i{i:.0f},j{j:.0f})'
        category = 'Forecast-SpringFlux'
    elif oname == 'recession':
        # Recession forecasts
        if 'flip_time' in longname or 'fliptime' in longname:
            component = 'Recession (flip_time)'
        else:
            component = 'Recession'
        category = 'Forecast-Recession'
    elif oname == 'budget':
        # Budget component forecasts
        if 'awanui' in longname:
            component = 'Awanui Stream'
            category = 'Budget-SW'
        elif 'poukawa' in longname:
            component = 'Poukawa Stream'
            category = 'Budget-SW'
        elif 'spring' in longname:
            component = 'Spring GHB'
            category = 'Budget-GHB'
        elif 'confined' in longname:
            component = 'Confined GHB'
            category = 'Budget-GHB'
        elif 'recharge' in longname:
            component = 'Recharge'
            category = 'Budget-Inflow'
        elif 'inflow' in longname:
            component = 'Total Inflow'
            category = 'Budget-Balance'
        elif 'stoss' in longname:
            component = 'Storage Change'
            category = 'Budget-Storage'
        elif 'sw' in longname and 'usecol:sw' in longname:
            component = 'Total SW'
            category = 'Budget-SW'
        elif 'total' in longname:
            component = 'Total Flux'
            category = 'Budget-Balance'
        elif 'inout' in longname:
            component = 'In-Out'
            category = 'Budget-Balance'
        elif 'percentdiscrepancy' in longname:
            component = 'Mass Balance Error'
            category = 'Budget-Balance'
        else:
            component = 'Other Budget'
            category = 'Budget-Other'
    else:
        component = 'Other'
        category = 'Other'

    # Check if it's a constraint
    is_constraint = 'greater' in obgnme or 'less' in obgnme
    constraint_type = None
    if 'greater' in obgnme:
        constraint_type = 'greater_than'
    elif 'less' in obgnme:
        constraint_type = 'less_than'

    return category, subcategory, component, period, is_constraint, constraint_type


def analyze_forecast_uncertainty(data):
    """Analyze forecast uncertainty reduction.

    Args:
        data: Dictionary from load_forecast_data

    Returns:
        pd.DataFrame: Forecast uncertainty analysis results
    """
    # Categories and components to exclude
    exclude_categories = ['Budget-Balance']
    exclude_components = ['Total Inflow', 'In-Out', 'Mass Balance Error', 'Total Flux']

    results = []

    for forecast_name in data['forecasts']:
        # Get forecast metadata
        forecast_row = data['forecast_obs'].loc[forecast_name]
        category, subcategory, component, period, is_constraint, constraint_type = categorize_forecast(forecast_row)

        # Skip excluded categories and components
        if category in exclude_categories or component in exclude_components:
            continue

        # Get prior and posterior ensembles
        prior_vals = data['prior'][forecast_name].values
        post_vals = data['posterior'][forecast_name].values

        # Calculate statistics
        prior_mean = np.mean(prior_vals)
        prior_std = np.std(prior_vals)
        prior_cv = (prior_std / abs(prior_mean) * 100) if prior_mean != 0 else np.nan

        post_mean = np.mean(post_vals)
        post_std = np.std(post_vals)
        post_cv = (post_std / abs(post_mean) * 100) if post_mean != 0 else np.nan

        # Uncertainty reduction
        std_reduction = ((prior_std - post_std) / prior_std * 100) if prior_std != 0 else np.nan
        cv_reduction = ((prior_cv - post_cv) / prior_cv * 100) if prior_cv != 0 and not np.isnan(prior_cv) else np.nan

        # Mean change
        mean_change = post_mean - prior_mean
        mean_change_pct = (mean_change / prior_mean * 100) if prior_mean != 0 else np.nan

        # Percentiles
        prior_p05 = np.percentile(prior_vals, 5)
        prior_p95 = np.percentile(prior_vals, 95)
        post_p05 = np.percentile(post_vals, 5)
        post_p95 = np.percentile(post_vals, 95)

        # Truth value if available
        truth_val = forecast_row['obsval']

        # Check if truth is within ensemble ranges
        truth_in_prior = (truth_val >= prior_p05) and (truth_val <= prior_p95)
        truth_in_post = (truth_val >= post_p05) and (truth_val <= post_p95)

        results.append({
            'forecast': forecast_name,
            'category': category,
            'subcategory': subcategory,
            'component': component,
            'period': period,
            'kper': forecast_row['kper'],
            'is_constraint': is_constraint,
            'constraint_type': constraint_type,
            'truth': truth_val,
            'prior_mean': prior_mean,
            'prior_std': prior_std,
            'prior_cv': prior_cv,
            'prior_p05': prior_p05,
            'prior_p95': prior_p95,
            'post_mean': post_mean,
            'post_std': post_std,
            'post_cv': post_cv,
            'post_p05': post_p05,
            'post_p95': post_p95,
            'mean_change': mean_change,
            'mean_change_pct': mean_change_pct,
            'std_reduction_pct': std_reduction,
            'cv_reduction_pct': cv_reduction,
            'truth_in_prior': truth_in_prior,
            'truth_in_post': truth_in_post
        })

    return pd.DataFrame(results)


def load_forecast_evolution(data):
    """Load forecast evolution across all iterations.

    Args:
        data: Dictionary from load_forecast_data

    Returns:
        pd.DataFrame: Forecast evolution data
    """
    # Categories and components to exclude
    exclude_categories = ['Budget-Balance']
    exclude_components = ['Total Inflow', 'In-Out', 'Mass Balance Error', 'Total Flux']

    evolution_data = []

    # Create forecast column list
    forecast_cols = ['real_name'] + list(data['forecasts'])

    for iter_num in data['iterations']:
        # Load observation ensemble for this iteration - only forecast columns for speed
        obs_file = os.path.join(data['ies_dir'], f"{data['model_name']}.{iter_num}.obs.csv")

        if not os.path.exists(obs_file):
            continue

        obs_df = pd.read_csv(obs_file, index_col=0, usecols=forecast_cols)
        forecast_vals = obs_df

        # Calculate statistics for each forecast
        for forecast_name in data['forecasts']:
            forecast_row = data['forecast_obs'].loc[forecast_name]
            category, subcategory, component, period, is_constraint, constraint_type = categorize_forecast(forecast_row)

            # Skip excluded categories and components
            if category in exclude_categories or component in exclude_components:
                continue

            vals = forecast_vals[forecast_name].values

            mean_val = np.mean(vals)
            std_val = np.std(vals)
            cv_val = (std_val / abs(mean_val) * 100) if mean_val != 0 else np.nan

            p05 = np.percentile(vals, 5)
            p50 = np.percentile(vals, 50)
            p95 = np.percentile(vals, 95)

            evolution_data.append({
                'iteration': iter_num,
                'forecast': forecast_name,
                'category': category,
                'component': component,
                'period': period,
                'kper': forecast_row['kper'],
                'mean': mean_val,
                'std': std_val,
                'cv': cv_val,
                'p05': p05,
                'p50': p50,
                'p95': p95,
                'truth': forecast_row['obsval']
            })

    return pd.DataFrame(evolution_data)


def create_visualizations(data, analysis_df, evolution_df, output_file):
    """Create comprehensive forecast analysis visualizations.

    Args:
        data: Dictionary from load_forecast_data
        analysis_df: DataFrame from analyze_forecast_uncertainty
        evolution_df: DataFrame from load_forecast_evolution
        output_file: Path to output figure file
    """
    fig = plt.figure(figsize=(20, 14))
    gs = gridspec.GridSpec(4, 3, figure=fig, hspace=0.35, wspace=0.3,
                          top=0.96, bottom=0.04, left=0.05, right=0.98)

    # Add title
    fig.suptitle(f'Forecast Analysis: {data["model_name"]}',
                 fontsize=16, fontweight='bold', y=0.99)

    # Color scheme
    category_colors = {
        'Surface Water': '#3498db',
        'GHB Flux': '#e74c3c',
        'Inflow': '#2ecc71',
        'Storage': '#f39c12',
        'Water Balance': '#9b59b6',
        'Other': '#95a5a6'
    }

    # 1. Uncertainty Reduction by Category
    ax1 = fig.add_subplot(gs[0, 0])
    category_stats = analysis_df.groupby('category').agg({
        'std_reduction_pct': 'mean',
        'cv_reduction_pct': 'mean'
    }).sort_values('std_reduction_pct', ascending=False)

    x = np.arange(len(category_stats))
    width = 0.35

    ax1.bar(x - width/2, category_stats['std_reduction_pct'], width,
            label='Std Dev Reduction', color='#3498db', alpha=0.7)
    ax1.bar(x + width/2, category_stats['cv_reduction_pct'], width,
            label='CV Reduction', color='#e74c3c', alpha=0.7)

    ax1.set_xlabel('Category', fontweight='bold')
    ax1.set_ylabel('Reduction (%)', fontweight='bold')
    ax1.set_title('Uncertainty Reduction by Category', fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(category_stats.index, rotation=45, ha='right')
    ax1.legend(frameon=False)
    ax1.grid(axis='y', alpha=0.3)
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8)

    # 2. Prior vs Posterior Uncertainty
    ax2 = fig.add_subplot(gs[0, 1])

    # Plot by component
    components = analysis_df.sort_values('prior_std', ascending=False)['component'].unique()[:8]
    comp_data = analysis_df[analysis_df['component'].isin(components)]

    x = np.arange(len(components))
    width = 0.35

    prior_stds = [comp_data[comp_data['component'] == c]['prior_std'].mean() for c in components]
    post_stds = [comp_data[comp_data['component'] == c]['post_std'].mean() for c in components]

    ax2.bar(x - width/2, prior_stds, width, label='Prior', color='#95a5a6', alpha=0.7)
    ax2.bar(x + width/2, post_stds, width, label='Posterior', color='#27ae60', alpha=0.7)

    ax2.set_xlabel('Component', fontweight='bold')
    ax2.set_ylabel('Standard Deviation', fontweight='bold')
    ax2.set_title('Prior vs Posterior Uncertainty (Top 8)', fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(components, rotation=45, ha='right', fontsize=9)
    ax2.legend(frameon=False)
    ax2.grid(axis='y', alpha=0.3)
    ax2.set_yscale('symlog')

    # 3. Mean Change
    ax3 = fig.add_subplot(gs[0, 2])

    # Plot mean change for top components
    comp_change = analysis_df.sort_values('mean_change', key=abs, ascending=False).head(10)

    colors = [category_colors.get(cat, '#95a5a6') for cat in comp_change['category']]

    ax3.barh(range(len(comp_change)), comp_change['mean_change'], color=colors, alpha=0.7)
    ax3.set_yticks(range(len(comp_change)))
    ax3.set_yticklabels([f"{c} ({p})" for c, p in zip(comp_change['component'], comp_change['period'])],
                        fontsize=9)
    ax3.set_xlabel('Mean Change', fontweight='bold')
    ax3.set_title('Mean Forecast Change (Top 10)', fontweight='bold')
    ax3.grid(axis='x', alpha=0.3)
    ax3.axvline(x=0, color='black', linestyle='-', linewidth=0.8)

    # 4. Forecast Evolution - Surface Water
    ax4 = fig.add_subplot(gs[1, :])

    sw_forecasts = evolution_df[evolution_df['category'] == 'Surface Water']

    for component in sw_forecasts['component'].unique():
        comp_data = sw_forecasts[sw_forecasts['component'] == component]

        for kper in sorted(comp_data['kper'].unique()):
            kper_data = comp_data[comp_data['kper'] == kper]

            label = f'{component} (kper{kper})'
            ax4.plot(kper_data['iteration'], kper_data['mean'], 'o-',
                    label=label, markersize=4, linewidth=1.5, alpha=0.8)
            ax4.fill_between(kper_data['iteration'], kper_data['p05'], kper_data['p95'],
                            alpha=0.2)

    # Add reinflation markers
    for r_iter in data['reinflate_iters']:
        ax4.axvline(x=r_iter, color='red', linestyle='--', linewidth=1.5, alpha=0.5)

    ax4.set_xlabel('Iteration', fontweight='bold')
    ax4.set_ylabel('Flux', fontweight='bold')
    ax4.set_title('Surface Water Forecast Evolution', fontweight='bold')
    ax4.legend(ncol=3, frameon=False, fontsize=8)
    ax4.grid(alpha=0.3)

    # 5. Forecast Evolution - GHB Flux
    ax5 = fig.add_subplot(gs[2, :])

    ghb_forecasts = evolution_df[evolution_df['category'] == 'GHB Flux']

    for component in ghb_forecasts['component'].unique():
        comp_data = ghb_forecasts[ghb_forecasts['component'] == component]

        for kper in sorted(comp_data['kper'].unique()):
            kper_data = comp_data[comp_data['kper'] == kper]

            label = f'{component} (kper{kper})'
            ax5.plot(kper_data['iteration'], kper_data['mean'], 's-',
                    label=label, markersize=4, linewidth=1.5, alpha=0.8)
            ax5.fill_between(kper_data['iteration'], kper_data['p05'], kper_data['p95'],
                            alpha=0.2)

    # Add reinflation markers
    for r_iter in data['reinflate_iters']:
        ax5.axvline(x=r_iter, color='red', linestyle='--', linewidth=1.5, alpha=0.5)

    ax5.set_xlabel('Iteration', fontweight='bold')
    ax5.set_ylabel('Flux', fontweight='bold')
    ax5.set_title('GHB Flux Forecast Evolution', fontweight='bold')
    ax5.legend(ncol=3, frameon=False, fontsize=8)
    ax5.grid(alpha=0.3)

    # 6. Water Balance Forecasts
    ax6 = fig.add_subplot(gs[3, 0])

    wb_forecasts = evolution_df[evolution_df['category'] == 'Water Balance']

    for component in wb_forecasts['component'].unique():
        comp_data = wb_forecasts[wb_forecasts['component'] == component]

        # Average across kpers
        comp_mean = comp_data.groupby('iteration')['mean'].mean()

        ax6.plot(comp_mean.index, comp_mean.values, 'o-',
                label=component, markersize=4, linewidth=1.5, alpha=0.8)

    # Add reinflation markers
    for r_iter in data['reinflate_iters']:
        ax6.axvline(x=r_iter, color='red', linestyle='--', linewidth=1.5, alpha=0.5)

    ax6.set_xlabel('Iteration', fontweight='bold')
    ax6.set_ylabel('Value', fontweight='bold')
    ax6.set_title('Water Balance Forecast Evolution', fontweight='bold')
    ax6.legend(frameon=False, fontsize=8)
    ax6.grid(alpha=0.3)
    ax6.set_yscale('symlog')

    # 7. Uncertainty Evolution (CV)
    ax7 = fig.add_subplot(gs[3, 1])

    # Show CV evolution for top 5 most uncertain forecasts
    final_iter = evolution_df['iteration'].max()
    final_cv = evolution_df[evolution_df['iteration'] == final_iter].sort_values('cv', ascending=False)
    top_forecasts = final_cv.head(5)['forecast'].tolist()

    for forecast_name in top_forecasts:
        fc_data = evolution_df[evolution_df['forecast'] == forecast_name]
        fc_row = data['forecast_obs'].loc[forecast_name]
        _, _, component, period, _, _ = categorize_forecast(fc_row)

        label = f'{component} ({period})'
        ax7.plot(fc_data['iteration'], fc_data['cv'], 'o-',
                label=label, markersize=4, linewidth=1.5, alpha=0.8)

    # Add reinflation markers
    for r_iter in data['reinflate_iters']:
        ax7.axvline(x=r_iter, color='red', linestyle='--', linewidth=1.5, alpha=0.5)

    ax7.set_xlabel('Iteration', fontweight='bold')
    ax7.set_ylabel('CV (%)', fontweight='bold')
    ax7.set_title('Uncertainty Evolution (Top 5)', fontweight='bold')
    ax7.legend(frameon=False, fontsize=8)
    ax7.grid(alpha=0.3)

    # 8. Truth vs Ensemble
    ax8 = fig.add_subplot(gs[3, 2])

    # Plot truth vs posterior mean
    ax8.scatter(analysis_df['truth'], analysis_df['post_mean'],
               c=[category_colors.get(cat, '#95a5a6') for cat in analysis_df['category']],
               alpha=0.6, s=60, edgecolors='black', linewidth=0.5)

    # Add 1:1 line
    lims = [
        np.min([ax8.get_xlim(), ax8.get_ylim()]),
        np.max([ax8.get_xlim(), ax8.get_ylim()]),
    ]
    ax8.plot(lims, lims, 'k-', alpha=0.5, zorder=0, linewidth=1)

    ax8.set_xlabel('Truth Value', fontweight='bold')
    ax8.set_ylabel('Posterior Mean', fontweight='bold')
    ax8.set_title('Truth vs Posterior Forecast', fontweight='bold')
    ax8.grid(alpha=0.3)
    ax8.set_aspect('equal')

    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Figure saved to: {output_file}")
    plt.close()


def save_summary_report(data, analysis_df, output_file):
    """Save markdown summary report.

    Args:
        data: Dictionary from load_forecast_data
        analysis_df: DataFrame from analyze_forecast_uncertainty
        output_file: Path to output markdown file
    """
    with open(output_file, 'w') as f:
        f.write(f"# Forecast Analysis Summary: {data['model_name']}\n\n")

        f.write("## Overview\n\n")
        f.write(f"Total forecasts: {len(data['forecasts'])}\n\n")
        f.write(f"Iterations analyzed: {len(data['iterations'])}\n\n")
        f.write(f"Last iteration: {data['last_iter']}\n\n")

        if data['reinflate_iters']:
            f.write(f"Reinflation iterations: {', '.join(map(str, data['reinflate_iters']))}\n\n")

        # Category summary
        f.write("## Forecast Categories\n\n")
        cat_summary = analysis_df.groupby('category').agg({
            'forecast': 'count',
            'std_reduction_pct': 'mean',
            'cv_reduction_pct': 'mean',
            'truth_in_post': 'sum'
        }).round(2)
        cat_summary.columns = ['Count', 'Avg Std Reduction %', 'Avg CV Reduction %', 'Truth in Post Range']
        f.write(cat_summary.to_markdown())
        f.write("\n\n")

        # Component summary
        f.write("## Top Uncertainty Reductions\n\n")
        top_reductions = analysis_df.sort_values('std_reduction_pct', ascending=False).head(10)
        summary_cols = ['component', 'period', 'prior_std', 'post_std', 'std_reduction_pct', 'truth_in_post']
        f.write(top_reductions[summary_cols].to_markdown(index=False))
        f.write("\n\n")

        # Largest mean changes
        f.write("## Largest Forecast Mean Changes\n\n")
        top_changes = analysis_df.iloc[analysis_df['mean_change'].abs().argsort()[::-1]].head(10)
        change_cols = ['component', 'period', 'prior_mean', 'post_mean', 'mean_change', 'mean_change_pct']
        f.write(top_changes[change_cols].to_markdown(index=False))
        f.write("\n\n")

        # Constraint forecasts
        constraints = analysis_df[analysis_df['is_constraint']]
        if len(constraints) > 0:
            f.write("## Constraint Forecasts\n\n")
            constraint_cols = ['component', 'period', 'constraint_type', 'truth', 'post_mean', 'post_std']
            f.write(constraints[constraint_cols].to_markdown(index=False))
            f.write("\n\n")

    print(f"Markdown summary saved to: {output_file}")


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description='Analyze PEST++ IES forecasts')
    parser.add_argument('model_name', nargs='?', default=None,
                       help='Model name (default: from setup.py)')
    parser.add_argument('--output-dir', default=None,
                       help='Output directory (default: models/{model_name})')

    args = parser.parse_args()

    # Get model name
    if args.model_name:
        model_name = args.model_name
    else:
        model_name = MODEL_NAME
        print(f"No model name provided, using MODEL_NAME from setup.py: {model_name}")

    # Set up paths
    model_dir = os.path.join('models', model_name)

    if not os.path.exists(model_dir):
        print(f"ERROR: Model directory not found: {model_dir}")
        return 1

    # Output directory
    if args.output_dir:
        output_dir = args.output_dir
    else:
        output_dir = model_dir

    os.makedirs(output_dir, exist_ok=True)

    print("="*80)
    print(f"PEST++ IES FORECAST ANALYSIS: {model_name}")
    print("="*80)
    print()

    # Load data
    print("Loading forecast data...")
    data = load_forecast_data(model_dir, model_name)

    print(f"Forecasts found: {len(data['forecasts'])}")
    print(f"Iterations: {len(data['iterations'])}")
    print(f"Reinflation iterations: {data['reinflate_iters']}")
    print()

    # Analyze uncertainty
    print("Analyzing forecast uncertainty...")
    analysis_df = analyze_forecast_uncertainty(data)

    # Load evolution
    print("Loading forecast evolution across iterations...")
    evolution_df = load_forecast_evolution(data)

    # Create visualizations
    print("Creating visualizations...")
    output_fig = os.path.join(output_dir, f'{model_name}_forecast_analysis.png')
    create_visualizations(data, analysis_df, evolution_df, output_fig)

    # Save results
    print("Saving results...")
    output_csv = os.path.join(output_dir, f'{model_name}_forecast_summary.csv')
    analysis_df.to_csv(output_csv, index=False)
    print(f"Forecast summary saved to: {output_csv}")

    output_evolution = os.path.join(output_dir, f'{model_name}_forecast_evolution.csv')
    evolution_df.to_csv(output_evolution, index=False)
    print(f"Forecast evolution saved to: {output_evolution}")

    output_md = os.path.join(output_dir, f'{model_name}_FORECAST_ANALYSIS.md')
    save_summary_report(data, analysis_df, output_md)

    print()
    print("="*80)
    print("ANALYSIS COMPLETE!")
    print("="*80)
    print()

    return 0


if __name__ == '__main__':
    sys.exit(main())
