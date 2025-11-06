import pandas as pd
import numpy as np

# Test to understand the data structure
prior_file = 'models/local_run24/pest/master_ies/local_run24.0.obs.csv'

# Read a small subset
df = pd.read_csv(prior_file, index_col=0, nrows=5)

print("DataFrame shape:", df.shape)
print("\nFirst few rows and columns:")
print(df.iloc[:5, :3])

print("\nColumn names (first 5):")
print(df.columns[:5].tolist())

print("\nIndex (realization names, first 5):")
print(df.index[:5].tolist())

# Check a specific observation column
ts_cols = [c for c in df.columns if 'ts-heads' in c and 'pk4_time' in c][:5]
print(f"\nFound {len(ts_cols)} ts-heads observations")
print("First 5 ts-heads columns:")
for col in ts_cols:
    print(f"  {col}")
    print(f"    Values: {df[col].values}")
