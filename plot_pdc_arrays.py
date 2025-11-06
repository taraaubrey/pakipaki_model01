"""
Plot prior data conflict arrays showing observation vs simulation conflicts
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re

# Read the PDC file
pdc_file = r"C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\local_run24\pest\master_ies_prior\local_run24.0.pdc.csv"
df = pd.read_csv(pdc_file)

# Parse the name column to extract i, j, kper, kstp
def parse_name(name_str):
    """Extract i, j, kper, kstp from name string"""
    # Pattern: oname:h-lyr1-kper{X}-kstp{Y}_otype:arr_i:{i}_j:{j}_zone:1
    kper_match = re.search(r'kper(\d+)', name_str)
    kstp_match = re.search(r'kstp(\d+)', name_str)
    i_match = re.search(r'_i:(\d+)', name_str)
    j_match = re.search(r'_j:(\d+)', name_str)

    return {
        'kper': int(kper_match.group(1)) if kper_match else None,
        'kstp': int(kstp_match.group(1)) if kstp_match else None,
        'i': int(i_match.group(1)) if i_match else None,
        'j': int(j_match.group(1)) if j_match else None
    }

# Parse all names
parsed_data = df['name'].apply(parse_name)
df_parsed = pd.DataFrame(parsed_data.tolist())
df = pd.concat([df, df_parsed], axis=1)

# Calculate the difference (obs_min - sim_min)
# Negative means sim < obs (simulation underestimates)
df['diff_obs_sim'] = df['obs_min'] - df['sim_min']

# Get unique kper/kstp combinations
unique_combos = df[['kper', 'kstp']].drop_duplicates().sort_values(['kper', 'kstp'])
print(f"Found {len(unique_combos)} unique kper/kstp combinations:")
print(unique_combos)

# Get grid dimensions
max_i = df['i'].max()
max_j = df['j'].max()
min_i = df['i'].min()
min_j = df['j'].min()
print(f"\nGrid extent: i=[{min_i}, {max_i}], j=[{min_j}, {max_j}]")

# Create output directory
output_dir = Path(r"C:\Users\tfo46\e_Python\e_projects\pakipaki_01\models\local_run24\figures\pdc_arrays")
output_dir.mkdir(parents=True, exist_ok=True)

# Plot each unique combination
for idx, (_, row) in enumerate(unique_combos.iterrows()):
    kper = int(row['kper'])
    kstp = int(row['kstp'])

    # Filter data for this combination
    mask = (df['kper'] == kper) & (df['kstp'] == kstp)
    subset = df[mask].copy()

    print(f"\nProcessing kper={kper}, kstp={kstp}: {len(subset)} observations")

    # Create arrays for plotting
    # Initialize with NaN
    diff_array = np.full((max_i + 1, max_j + 1), np.nan)
    sim_min_array = np.full((max_i + 1, max_j + 1), np.nan)
    sim_max_array = np.full((max_i + 1, max_j + 1), np.nan)

    # Fill arrays with data
    for _, data_row in subset.iterrows():
        i = int(data_row['i'])
        j = int(data_row['j'])
        diff_array[i, j] = data_row['diff_obs_sim']
        sim_min_array[i, j] = data_row['sim_min']
        sim_max_array[i, j] = data_row['sim_max']

    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f'Prior Data Conflict - kper={kper}, kstp={kstp}', fontsize=14, fontweight='bold')

    # Plot 1: Difference (obs_min - sim_min)
    im1 = axes[0].imshow(diff_array, cmap='RdBu_r', origin='upper', aspect='auto')
    axes[0].set_title('obs_min - sim_min\n(negative = sim underestimates)')
    axes[0].set_xlabel('j index')
    axes[0].set_ylabel('i index')
    cbar1 = plt.colorbar(im1, ax=axes[0])
    cbar1.set_label('Head difference (m)', rotation=270, labelpad=20)

    # Plot 2: sim_min
    im2 = axes[1].imshow(sim_min_array, cmap='viridis', origin='upper', aspect='auto')
    axes[1].set_title('sim_min\n(minimum simulated head)')
    axes[1].set_xlabel('j index')
    axes[1].set_ylabel('i index')
    cbar2 = plt.colorbar(im2, ax=axes[1])
    cbar2.set_label('Head (m)', rotation=270, labelpad=20)

    # Plot 3: sim_max
    im3 = axes[2].imshow(sim_max_array, cmap='viridis', origin='upper', aspect='auto')
    axes[2].set_title('sim_max\n(maximum simulated head)')
    axes[2].set_xlabel('j index')
    axes[2].set_ylabel('i index')
    cbar3 = plt.colorbar(im3, ax=axes[2])
    cbar3.set_label('Head (m)', rotation=270, labelpad=20)

    plt.tight_layout()

    # Save figure
    output_file = output_dir / f'pdc_kper{kper}_kstp{kstp}.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_file}")

    plt.close()

print(f"\nAll plots saved to: {output_dir}")

# Create a summary figure showing statistics
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Prior Data Conflict Summary Statistics', fontsize=14, fontweight='bold')

# Group by kper/kstp for summary stats
summary = df.groupby(['kper', 'kstp']).agg({
    'diff_obs_sim': ['mean', 'std', 'min', 'max'],
    'distance': ['mean', 'max']
}).reset_index()
summary.columns = ['_'.join(col).strip('_') for col in summary.columns.values]

# Plot 1: Mean difference by time
ax = axes[0, 0]
summary['time_label'] = summary.apply(lambda x: f"k{int(x['kper'])}-s{int(x['kstp'])}", axis=1)
ax.bar(range(len(summary)), summary['diff_obs_sim_mean'], yerr=summary['diff_obs_sim_std'])
ax.set_xticks(range(len(summary)))
ax.set_xticklabels(summary['time_label'], rotation=45, ha='right')
ax.set_ylabel('Mean obs-sim difference (m)')
ax.set_title('Mean bias over time')
ax.axhline(0, color='k', linestyle='--', alpha=0.3)
ax.grid(True, alpha=0.3)

# Plot 2: Range of differences
ax = axes[0, 1]
ax.plot(range(len(summary)), summary['diff_obs_sim_min'], 'b-o', label='Min difference')
ax.plot(range(len(summary)), summary['diff_obs_sim_max'], 'r-o', label='Max difference')
ax.set_xticks(range(len(summary)))
ax.set_xticklabels(summary['time_label'], rotation=45, ha='right')
ax.set_ylabel('Difference (m)')
ax.set_title('Range of obs-sim differences')
ax.axhline(0, color='k', linestyle='--', alpha=0.3)
ax.legend()
ax.grid(True, alpha=0.3)

# Plot 3: Distance metric
ax = axes[1, 0]
ax.bar(range(len(summary)), summary['distance_mean'])
ax.set_xticks(range(len(summary)))
ax.set_xticklabels(summary['time_label'], rotation=45, ha='right')
ax.set_ylabel('Mean distance')
ax.set_title('Mean conflict distance')
ax.grid(True, alpha=0.3)

# Plot 4: Number of conflicts
ax = axes[1, 1]
conflict_counts = df.groupby(['kper', 'kstp']).size().reset_index(name='count')
conflict_counts['time_label'] = conflict_counts.apply(lambda x: f"k{int(x['kper'])}-s{int(x['kstp'])}", axis=1)
ax.bar(range(len(conflict_counts)), conflict_counts['count'])
ax.set_xticks(range(len(conflict_counts)))
ax.set_xticklabels(conflict_counts['time_label'], rotation=45, ha='right')
ax.set_ylabel('Number of conflicts')
ax.set_title('Conflict count by time')
ax.grid(True, alpha=0.3)

plt.tight_layout()
summary_file = output_dir / 'pdc_summary.png'
plt.savefig(summary_file, dpi=150, bbox_inches='tight')
print(f"\nSummary plot saved: {summary_file}")
plt.close()

# Print summary statistics
print("\n" + "="*60)
print("SUMMARY STATISTICS")
print("="*60)
print(f"\nTotal observations with conflicts: {len(df)}")
print(f"\nOverall statistics:")
print(f"  Mean obs-sim difference: {df['diff_obs_sim'].mean():.3f} m")
print(f"  Std obs-sim difference:  {df['diff_obs_sim'].std():.3f} m")
print(f"  Min obs-sim difference:  {df['diff_obs_sim'].min():.3f} m (sim >> obs)")
print(f"  Max obs-sim difference:  {df['diff_obs_sim'].max():.3f} m (obs >> sim)")
print(f"\nNegative differences (sim < obs): {(df['diff_obs_sim'] < 0).sum()} / {len(df)} ({100*(df['diff_obs_sim'] < 0).sum()/len(df):.1f}%)")
print(f"Positive differences (sim > obs): {(df['diff_obs_sim'] > 0).sum()} / {len(df)} ({100*(df['diff_obs_sim'] > 0).sum()/len(df):.1f}%)")

print("\n" + "="*60)
