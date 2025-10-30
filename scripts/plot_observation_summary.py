"""
Standalone script to generate observation summary table
Similar to parameter summary but for observations
"""
import os
import pandas as pd

# Configuration
run_name = 'local_run3'
m_d = f'models/{run_name}/pest/master_ies'
output = f'models/{run_name}/figures'
obs_data_fn = os.path.join(m_d, f"{run_name}.obs_data.csv")

# Read observation data
obs_data = pd.read_csv(obs_data_fn)

print(f"\nLoaded {len(obs_data)} observations from {obs_data_fn}")
print(f"\nObservation groups found: {obs_data['obgnme'].unique().tolist()}")

# Create observation summary
obs_summary_rows = []

# Group observations and create summary
for obgnme, group in obs_data.groupby('obgnme'):
    print(f"\nProcessing group: {obgnme} ({len(group)} observations)")

    # Count observations with non-zero weights
    count_weighted = (group['weight'] > 0).sum()
    count_total = len(group)

    # Get standard deviation (mean of non-zero values)
    std = group.loc[group['standard_deviation'] > 0, 'standard_deviation'].mean()

    # Get adjusted weight (mean of non-zero values)
    adj_weight = group.loc[group['weight'] > 0, 'weight'].mean()

    # Determine observation type and details
    if obgnme == 'ts-heads':
        # Get unique columns
        unique_cols = group['usecol'].unique()
        for col in unique_cols:
            col_group = group[group['usecol'] == col]
            count = (col_group['weight'] > 0).sum()
            std_col = col_group.loc[col_group['standard_deviation'] > 0, 'standard_deviation'].mean()
            weight_col = col_group.loc[col_group['weight'] > 0, 'weight'].mean()

            # Get time range
            times = col_group.loc[col_group['weight'] > 0, 'time']
            if len(times) > 0:
                time_range = f"{int(times.min())}-{int(times.max())}"
            else:
                time_range = "-"

            # Determine name and type
            if col == 'pk4':
                name = 'PK4'
                obs_type = 'head'
            elif col == 'pk4-spr-diff':
                name = 'PK4 - Spring'
                obs_type = 'head difference'
            elif col == 'pk4-conf-diff':
                name = 'PK4 - Confined'
                obs_type = 'head difference'
            else:
                name = col
                obs_type = 'head'

            obs_summary_rows.append({
                'Group': obgnme,
                'Name': name,
                'Type': obs_type,
                'Stress period': time_range,
                'Inequality constraint': '-',
                'Count': count,
                'Std': f"{std_col:.3g}" if pd.notna(std_col) else '-',
                'Adjusted weight': f"{weight_col:.3g}" if pd.notna(weight_col) else '-'
            })

    elif obgnme == 'AWq':
        # Get kper range for non-zero weights
        kpers = group.loc[group['weight'] > 0, 'kper']
        if len(kpers) > 0:
            kper_range = f"{int(kpers.min())}"
        else:
            kper_range = "-"

        obs_summary_rows.append({
            'Group': obgnme,
            'Name': 'Awanui cell flux',
            'Type': 'flux',
            'Stress period': kper_range,
            'Inequality constraint': '-',
            'Count': count_weighted,
            'Std': f"{std:.3g}" if pd.notna(std) else '-',
            'Adjusted weight': f"{adj_weight:.3g}" if pd.notna(adj_weight) else '-'
        })

    elif obgnme == 'recession':
        # Get unique columns
        unique_cols = group['usecol'].unique()
        for col in unique_cols:
            col_group = group[group['usecol'] == col]
            count = (col_group['weight'] > 0).sum()
            std_col = col_group.loc[col_group['standard_deviation'] > 0, 'standard_deviation'].mean()
            weight_col = col_group.loc[col_group['weight'] > 0, 'weight'].mean()

            # Determine name
            if col == 'pk4':
                name = 'PK4'
            elif col == 'pk4-spr-diff':
                name = 'PK4 - Spring'
            elif col == 'pk4-conf-diff':
                name = 'PK4 - Confined'
            else:
                name = col

            obs_summary_rows.append({
                'Group': obgnme,
                'Name': name,
                'Type': 'recession rate',
                'Stress period': '-',
                'Inequality constraint': '-',
                'Count': count,
                'Std': f"{std_col:.3g}" if pd.notna(std_col) else '-',
                'Adjusted weight': f"{weight_col:.3g}" if pd.notna(weight_col) else '-'
            })

    elif obgnme.startswith('less_') or obgnme.startswith('greater_'):
        # Budget or head-arr with constraints
        base_group = obgnme.replace('less_', '').replace('greater_', '')
        constraint_type = '< 0' if obgnme.startswith('less_') else '> 0'

        if base_group == 'budget':
            # Get unique columns
            unique_cols = group['usecol'].unique()
            for col in unique_cols:
                col_group = group[group['usecol'] == col]
                count = (col_group['weight'] > 0).sum()
                std_col = col_group.loc[col_group['standard_deviation'] > 0, 'standard_deviation'].mean()
                weight_col = col_group.loc[col_group['weight'] > 0, 'weight'].mean()

                # Get kper range
                kpers = col_group.loc[col_group['weight'] > 0, 'kper']
                if len(kpers) > 0:
                    kper_range = f"{int(kpers.min())}-{int(kpers.max())}"
                else:
                    kper_range = "-"

                obs_summary_rows.append({
                    'Group': 'budget',
                    'Name': f"{col.capitalize()} Flux",
                    'Type': 'flux',
                    'Stress period': kper_range,
                    'Inequality constraint': constraint_type,
                    'Count': count,
                    'Std': f"{std_col:.3g}" if pd.notna(std_col) else '-',
                    'Adjusted weight': f"{weight_col:.3g}" if pd.notna(weight_col) else '-'
                })
        elif base_group == 'head-arr':
            # Get kper range
            # Note: head-arr observations don't have kper directly, need to parse from oname
            kpers_set = set()
            for oname in group['oname'].unique():
                if 'kper' in oname:
                    kper_str = oname.split('kper')[1].split('-')[0].split('_')[0]
                    try:
                        kpers_set.add(int(kper_str))
                    except:
                        pass

            if kpers_set:
                kper_min, kper_max = min(kpers_set), max(kpers_set)
                if kper_min == kper_max:
                    kper_range = f"{kper_min}"
                else:
                    kper_range = f"{kper_min}-{kper_max}"
            else:
                kper_range = "-"

            obs_summary_rows.append({
                'Group': 'head-arr',
                'Name': 'Heads (array)',
                'Type': 'head',
                'Stress period': kper_range,
                'Inequality constraint': '< top + 1m',
                'Count': count_weighted,
                'Std': f"{std:.3g}" if pd.notna(std) else '-',
                'Adjusted weight': f"{adj_weight:.3g}" if pd.notna(adj_weight) else '-'
            })

# Create DataFrame
obs_summary = pd.DataFrame(obs_summary_rows)

# Save to CSV
output_file = os.path.join(output, 'output.observation_summary.csv')
obs_summary.to_csv(output_file, index=False)

print(f"\n{'='*80}")
print(f"Observation Summary saved to: {output_file}")
print(f"{'='*80}\n")
print(obs_summary.to_string(index=False))
print(f"\n{'='*80}")
print("Done!")
