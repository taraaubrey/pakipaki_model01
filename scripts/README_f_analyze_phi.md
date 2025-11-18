# f_analyze_phi.py - PEST++ IES Phi Analysis Tool

## Purpose

Analyzes PEST++ IES phi results to understand observation group performance, particularly focusing on arr-h (array-based head) observations. Run this after completing `e_pest_IES_HM.py`.

## Usage

```bash
# Use the MODEL_NAME from setup.py (default)
python scripts/f_analyze_phi.py

# Or specify a specific model run
python scripts/f_analyze_phi.py local_run33
```

## What It Does

1. **Reads phi results** from `models/{MODEL_NAME}/pest/master_ies/`
2. **Identifies observation groups** by type:
   - **ARR-H**: Array-based head observations (arr-awq, arr-spq, arr-confq, arr-h)
   - **Head**: Point head observations (ts-heads)
   - **Recession**: Recession rate observations
   - **Flux**: Flux/flow observations
   - **Budget**: Budget observations
3. **Calculates performance metrics**:
   - Phi reduction percentage from initial to final
   - Final contribution to total residual phi
   - Phi evolution over iterations
4. **Detects reinflation iterations** (where phi increases)

## Outputs

### 1. Console Output
- Observation group identification with weights
- Performance summary table
- Final phi contributions by type
- Total phi evolution

### 2. PNG Figure: `{MODEL_NAME}_phi_analysis.png`
Comprehensive visualization with:
- Total phi evolution with reinflation markers
- Phi evolution by group type (ARR-H, Head, Recession)
- Normalized phi comparison
- Final phi values by group
- Phi reduction percentages

### 3. CSV Table: `{MODEL_NAME}_phi_summary.csv`
All observation groups with:
- Type, Group, Obs_Name, Weight
- Initial_Phi, Final_Phi, Reduction_%
- Final_Contribution_%

### 4. Markdown Summary: `{MODEL_NAME}_PHI_ANALYSIS.md`
- Overview with total iterations and phi reduction
- Reinflation iterations identified
- ARR-H performance table
- Complete observation group table

## Interpreting Results

### ARR-H Performance

**Good performance indicators:**
- Phi reduction >99%
- Final contribution to total phi <5%
- Steady downward trend (excluding reinflation)

**Potential issues:**
- Phi reduction <90%
- High final contribution (>20%)
- Persistent upward trend after reinflation recovery

### Reinflation (Iteration 4 in local_run33)

Reinflation is **normal** and expected in PEST++ IES:
- The algorithm intentionally increases phi temporarily
- Helps maintain ensemble variance
- Should recover in subsequent iterations
- Marked with red dashed lines in plots

### Final Phi Contributions

Shows which observation types dominate residual error:
- **Low (<2%)**: Well-matched observations
- **Medium (2-10%)**: Acceptable fit
- **High (>20%)**: Primary contributors to mismatch

For local_run33:
- ARR-H: 1.21% (excellent)
- Head: 0.08% (excellent)
- Recession: 96.48% ⚠️ (primary concern)
- Budget: 2.23% (acceptable)

## Example Interpretation (local_run33)

```
ARR-H observations achieved 99.6-99.98% phi reduction
Total ARR-H contribution to final phi: 1.21%

✓ ARR-H observations are performing well
✓ Reinflation at iteration 4 is normal and recovered by iteration 8
⚠ Recession observations (71-78% reduction, 96% of total phi) need attention
```

## Integration into Workflow

Add to your workflow after `e_pest_IES_HM.py`:

```python
# In your main workflow script
from setup import MODEL_NAME
import os

# ... run a_model_setup.py through e_pest_IES_HM.py ...

# Analyze phi results
print("Analyzing phi results...")
os.system(f"python scripts/f_analyze_phi.py {MODEL_NAME}")
```

Or run manually after each PEST++ IES run to check convergence and observation performance.
