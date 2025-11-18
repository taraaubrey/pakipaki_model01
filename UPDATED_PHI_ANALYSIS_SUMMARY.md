# Updated Phi Analysis Script - Summary

## Changes Made

### 1. **Added usecol to Observation Names**

Observation names now include the `usecol` field for better differentiation:

**Before:**
- `ts-heads` (ambiguous - which head type?)
- `recession` (ambiguous - which recession metric?)
- `budget` (ambiguous - which budget component?)

**After:**
- `ts-heads:pk4` (heads at Poukawa location 4)
- `ts-heads:pk4-conf-diff` (heads difference from confined layer)
- `ts-heads:pk4-spr-diff` (heads difference from spring layer)
- `recession:pk4` (recession rate at pk4)
- `recession:pk4-conf-diff` (recession rate for confined difference)
- `budget:confined` (budget for confined area)
- `budget:inflow` (budget inflow component)
- `arr-awq:awq` (array heads for Awanui)
- `arr-spq:spq` (array heads for springs)
- `arr-confq:confq` (array heads for confined area)

### 2. **Bar Graph Now Shows Only Non-Zero Weighted Observations**

The bottom bar graph (Plot 7: "Phi Reduction from Initial") now filters out observations with weight=0, making it cleaner and more focused on active observations.

**Result for local_run33:**
- Shows only 12 observations (instead of 25)
- All displayed observations are actively contributing to calibration
- Title updated to: "Phi Reduction from Initial by Observation Group (Non-Zero Weights Only)"

## Current Results for local_run33

### ARR-H Observations (Now with usecol)

| Group | Full Name | Weight | Reduction | Contribution |
|-------|-----------|--------|-----------|--------------|
| less_16obg22 | arr-h | 0.30 | 99.98% | 0.13% |
| 19obg5 | **arr-spq:spq** | 0.05 | 99.97% | 0.04% |
| 18obg4 | **arr-awq:awq** | 0.05 | 99.95% | 0.07% |
| greater_20obg6 | **arr-confq:confq** | 0.10 | 99.63% | 0.98% |

**Total ARR-H contribution: 1.21%** ✅

### Head Observations (Now Differentiated)

| Group | Full Name | Weight | Reduction | Contribution |
|-------|-----------|--------|-----------|--------------|
| 13obg2 | **ts-heads:pk4** | 0.30 | 99.9968% | 0.026% |
| 1obg0 | **ts-heads:pk4-conf-diff** | 0.30 | 99.9968% | 0.026% |
| 2obg1 | **ts-heads:pk4-spr-diff** | 0.30 | 99.9967% | 0.026% |

### Recession Observations (Now Differentiated)

| Group | Full Name | Weight | Reduction | Contribution |
|-------|-----------|--------|-----------|--------------|
| 3obg10 | **recession:pk4** | 0.05 | 78.30% | 29.13% |
| 22obg8 | **recession:pk4-conf-diff** | 0.05 | 78.30% | 29.37% |
| 23obg9 | **recession:pk4-spr-diff** | 0.05 | 71.94% | 37.98% |

**Total Recession contribution: 96.48%** ⚠️ *These are the primary source of residual phi*

### Budget Observations (Now Differentiated)

| Group | Full Name | Weight | Reduction | Contribution |
|-------|-----------|--------|-----------|--------------|
| greater_5obg12 | **budget:confined** | 0.20 | 99.83% | 0.94% |
| less_6obg13 | **budget:inflow** | 0.20 | 99.76% | 1.28% |

## Key Insights

### ARR-H Performance is Excellent
- All arr-h groups achieve >99.6% phi reduction
- Combined contribution to total residual phi: **1.21%** (minimal)
- The spike at iteration 4 (reinflation) is **normal and expected**

### The Real Issue: Recession Rates
All three recession observation groups struggle:
1. **recession:pk4-spr-diff** (23obg9): 71.94% reduction → 37.98% of total phi
2. **recession:pk4-conf-diff** (22obg8): 78.30% reduction → 29.37% of total phi
3. **recession:pk4** (3obg10): 78.30% reduction → 29.13% of total phi

**Together they account for 96.48% of all residual phi**

This suggests:
- Recession rate observations are difficult to match with current parameterization
- Parameters controlling transient head decline may need tighter constraints
- Storage properties (specific storage, specific yield) may need review
- Consider whether recession rate observation weights are appropriate

## Using the Updated Script

```bash
# Default (uses MODEL_NAME from setup.py)
python scripts/f_analyze_phi.py

# For any specific model run
python scripts/f_analyze_phi.py local_run34
```

## Outputs

All outputs now include the improved `usecol` differentiation:

1. **Console output** - Shows grouped observations with full names
2. **PNG figure** - `{MODEL_NAME}_phi_analysis.png` with bar graph showing only weighted obs
3. **CSV table** - `{MODEL_NAME}_phi_summary.csv` with full observation names
4. **Markdown report** - `{MODEL_NAME}_PHI_ANALYSIS.md` with detailed tables

## Recommendations

Based on the updated analysis:

1. ✅ **ARR-H observations are working well** - No changes needed

2. ⚠️ **Focus on recession observations:**
   - Review if recession rate calculation is correct
   - Check if observation weights (0.05) are appropriate
   - Consider if temporal weighting scheme is needed
   - Examine storage parameter bounds (STO package)
   - Review if recession observations are physically realistic

3. 📊 **Use this script after each PEST++ run** to track convergence and identify problem observations

## Script Location

- **Main script**: `scripts/f_analyze_phi.py`
- **Documentation**: `scripts/README_f_analyze_phi.md`
