# Stuck Realizations Diagnosis - local_run33

## Problem Summary

**39.3% of realizations (86 out of 219) show minimal improvement:**
- Less than 50% phi reduction
- Some actually got worse (negative reduction!)
- Worst case: Realization 226 got 1199% WORSE (phi increased 13×)

## Root Cause Analysis

### 🔍 Finding #1: Stuck Realizations Are Under-Updated

**Parameter Changes (Mean Absolute Change):**
- **Stuck realizations**: 20,029 (mean change)
- **Good realizations**: 58,950 (mean change)
- **Ratio**: 0.34× (stuck realizations change 3× LESS)

**Interpretation:**
✅ Stuck realizations ARE being updated (86% of parameters change)
⚠️ BUT they're being updated with **3× smaller magnitude** than good realizations

This suggests the PEST++ IES algorithm is:
- Giving smaller parameter updates to these realizations
- Possibly due to low ensemble variance
- Or these realizations have poor sensitivity

### 🔍 Finding #2: Parameters Are NOT Hitting Bounds

**Parameters at Bounds:**
- Stuck realizations: Only 1.30% at bounds
- This is very low - NOT the issue

**Interpretation:**
✅ Parameter bounds are NOT the problem
✅ Parameters have room to move

### 🔍 Finding #3: Worst Offenders

| Realization | Initial Phi | Final Phi | Change | Status |
|-------------|-------------|-----------|--------|--------|
| **226** | 3.44e+13 | **4.46e+14** | **+1199%** | 🔴 Got 13× WORSE |
| **7** | 2.27e+14 | **5.12e+14** | **+126%** | 🔴 Got 2× worse |
| **170** | 8.44e+13 | **1.28e+14** | **+51%** | 🔴 Got worse |
| **17** | 3.15e+13 | **4.72e+13** | **+50%** | 🔴 Got worse |
| **146** | 4.93e+14 | **6.94e+14** | **+41%** | 🔴 Got worse |

## Why This Is Happening

### Theory 1: Low Ensemble Variance (Most Likely) ⭐

**What's happening:**
- Initial ensemble doesn't have enough diversity
- All realizations cluster around similar parameter values
- PEST++ can't find good directions to update stuck realizations
- Updates are conservative (small magnitude)

**Evidence:**
- 86% of parameters change, but magnitude is 3× smaller
- Similar percentage of parameters changed (86% vs 84.5%)
- But absolute changes are much smaller

**Solution:**
```python
# In d_pest_IES_prior.py or e_pest_IES_HM.py, increase ensemble diversity:

# Option 1: Increase initial parameter variance
pst.parameter_data.loc[:, 'parval1'] = prior_mean
pst.parameter_data.loc[:, 'parstd'] = prior_std * 2.0  # Increase from 1.0 to 2.0

# Option 2: Use larger prior ensemble
NREALS_PRIOR = 500  # Increase from 250

# Option 3: Loosen parameter bounds to allow more exploration
```

### Theory 2: Poor Observation-Parameter Sensitivity

**What's happening:**
- Some realizations land in parameter space where observations are insensitive
- Changes don't improve phi because observations don't respond
- Algorithm gives up on updating these realizations strongly

**Evidence:**
- Some realizations got worse (bad parameter updates)
- Updates are small (algorithm is uncertain)

**Solution:**
```python
# In PEST++ control file, adjust localization:
ies_localize_how(empirical)  # Use empirical localization
ies_num_reals(100)  # Ensure enough realizations for sensitivity
```

### Theory 3: Reinflation Issues

**What's happening:**
- Iteration 4 is a reinflation iteration
- Some realizations get bad parameter combinations during reinflation
- Never recover

**Evidence:**
- All runs show phi spike at iteration 4
- Some realizations may have been "ruined" by reinflation

**Solution:**
```python
# In PEST++ control file:
ies_no_noise(true)  # Disable reinflation noise
# OR
ies_lambda_mults(0.1, 1.0, 10.0)  # More conservative lambda schedule
```

## Recommended Actions

### Immediate Actions (for next run):

1. **Increase ensemble diversity:**
   ```python
   # In scripts/d_pest_IES_prior.py
   # When generating prior ensemble, increase variance:
   prior_std_multiplier = 2.0  # Double the prior standard deviation
   ```

2. **Reduce recession observation weights:**
   ```python
   # As discussed earlier, reduce from 0.05 to 0.015
   # This may help realizations that are stuck on recession
   ```

3. **Try more conservative lambda:**
   ```python
   # In PEST++ control file:
   ies_lambda_mults(0.5, 1.0, 2.0)  # Less aggressive than default
   ```

### Medium-term Actions:

4. **Increase prior ensemble size:**
   ```python
   NREALS_PRIOR = 500  # Up from 250
   # More realizations = better ensemble variance
   ```

5. **Check parameter correlation structure:**
   ```python
   # Review geostatistical parameters in parameterization:
   # - Variogram range
   # - Correlation structure
   # May need more spatial variability
   ```

6. **Use subset selection:**
   ```python
   # In PEST++ control:
   ies_subset_size(100)  # Use best 100 realizations only
   # Drops poorly performing realizations
   ```

### Diagnostic Actions:

7. **Analyze ensemble variance evolution:**
   ```python
   # Check if ensemble variance is collapsing
   # Plot parameter variance across iterations
   # Should stay > 0.1 * initial variance
   ```

8. **Review localization:**
   ```python
   # Current localization may be too strong
   # Try looser localization or adaptive localization
   ```

## Quick Test

To test if low ensemble variance is the issue, try this:

```python
# In your next run, after generating prior ensemble:
# Check the variance of key parameters

import pandas as pd
par0 = pd.read_csv('models/{MODEL_NAME}/pest/master_ies/{MODEL_NAME}.0.par.csv')

# Check variance for a few key parameters
key_params = ['p0', 'p100', 'p500']  # Sample some parameters
for param in key_params:
    var = par0[param].var()
    mean = par0[param].mean()
    cv = var**0.5 / mean if mean != 0 else 0
    print(f"{param}: mean={mean:.3f}, std={var**0.5:.3f}, CV={cv:.3f}")

# If CV < 0.5, you have low variance problem
```

## Expected Improvement

If ensemble variance is the issue (most likely), you should see:

**Before fix:**
- 39% stuck realizations
- Mean parameter change ratio: 0.34× (stuck vs good)

**After fix (increased variance):**
- <10% stuck realizations
- Mean parameter change ratio: >0.7× (stuck vs good)
- More realizations achieving >90% phi reduction

## Files Generated

- `models/local_run33/stuck_realizations.csv` - List of 86 stuck realizations
- `models/local_run33/realization_performance.csv` - All realizations with phi reduction

## Summary

**Problem**: 39% of realizations barely improve (3× smaller parameter updates)

**Root Cause**: Low ensemble variance - realizations too similar

**Solution**: Increase prior parameter variance by 2× and consider:
- Larger ensemble size (500 vs 250)
- More conservative lambda schedule
- Reduced recession weights (0.015 vs 0.05)
- Subset selection to drop poor realizations

**Test**: Check if parameter CV < 0.5 in initial ensemble
