# Observation Weight Adjustment Recommendations

## Current Situation (local_run33)

### Contribution to Total Residual Phi

| Type | Current Weight | Final Contribution % | Problem |
|------|----------------|---------------------|---------|
| **Recession** | 0.05 | **96.48%** | Dominates visualization |
| ARR-H | 0.05-0.30 | 1.21% | Hard to see |
| Flux (budget:inflow) | 0.20 | 1.28% | Hard to see |
| Budget | 0.20 | 0.94% | Hard to see |
| Head | 0.30 | 0.08% | Barely visible |

## Recommendation Options

### Option 1: Moderate Reduction (Recommended) ⭐
**Goal**: Make other groups visible while keeping recession important

**New Weights:**
```python
# In phi_factors.csv or parameterization setup
recession weights: 0.05 → 0.015  # Reduce by ~70%
```

**Expected Result:**
- Recession contribution: ~30% (down from 96%)
- ARR-H contribution: ~4% (up from 1.2%)
- Makes trends visible without ignoring recession

**Pros:**
- ✅ Other observation groups become visible in plots
- ✅ Recession still important enough to drive calibration
- ✅ Balanced view of all observation types

**Cons:**
- ⚠️ Recession observations get less influence on parameters
- ⚠️ May converge slower on recession metrics

### Option 2: Aggressive Reduction
**Goal**: Prioritize well-performing observations

**New Weights:**
```python
recession weights: 0.05 → 0.005  # Reduce by 90%
```

**Expected Result:**
- Recession contribution: ~10% (down from 96%)
- ARR-H contribution: ~12% (up from 1.2%)
- Other observations dominate calibration

**Pros:**
- ✅ Excellent visualization of well-performing groups
- ✅ Parameters optimized for arr-h and heads
- ✅ May improve overall phi if recession is unattainable

**Cons:**
- ❌ Recession observations mostly ignored
- ❌ Transient behavior may not be captured well
- ❌ Risk of overfitting to head observations

### Option 3: Selective Adjustment
**Goal**: Keep some recession weights, reduce problematic ones

**New Weights:**
```python
# recession:pk4 (best performer) → keep at 0.05
# recession:pk4-conf-diff → reduce to 0.01
# recession:pk4-spr-diff (worst) → reduce to 0.005
```

**Expected Result:**
- Focus on most achievable recession metric
- Still attempts to match transient behavior
- Better balance across observation types

**Pros:**
- ✅ Nuanced approach targeting poor performers
- ✅ Keeps recession that's working better
- ✅ Moderate improvement in visualization

**Cons:**
- ⚠️ More complex to tune
- ⚠️ May still be recession-dominated

### Option 4: Two-Stage Approach (Advanced) 🎯
**Goal**: First calibrate to heads, then add recession

**Stage 1 Weights:**
```python
recession weights: 0.001  # Minimal
arr-h, heads, budget: Keep current (0.05-0.30)
```

**Stage 2 Weights (after heads converge):**
```python
recession weights: 0.02  # Moderate
arr-h, heads: Keep current
```

**Pros:**
- ✅ Best visualization during initial calibration
- ✅ Parameters first match spatial heads well
- ✅ Then refine for transient behavior
- ✅ Can run as separate PEST++ runs

**Cons:**
- ❌ Requires multiple PEST++ runs
- ❌ More workflow complexity

## Practical Implementation

### How to Change Weights

**Method 1: Edit phi_factors.csv directly**
```bash
# After running c_parameterization.py, before e_pest_IES_HM.py
# Edit: models/{MODEL_NAME}/pest/{MODEL_NAME}_template/phi_factors.csv

# Change lines like:
22obg8,0.05  → 22obg8,0.015
23obg9,0.05  → 23obg9,0.015
3obg10,0.05  → 3obg10,0.015
```

**Method 2: Modify in parameterization code**
Edit `scripts/c_parameterization.py` (or equivalent) where phi factors are set:

```python
# Find where phi_factors are defined
phi_factors = {
    '1obg0': 0.3,
    '2obg1': 0.3,
    '13obg2': 0.3,
    # ... other groups ...
    '3obg10': 0.015,   # Changed from 0.05
    '22obg8': 0.015,   # Changed from 0.05
    '23obg9': 0.015,   # Changed from 0.05
}
```

## What to Expect After Changing Weights

### Immediate Effects:
1. **Phi will increase initially** - Recession has less weight, so total phi reflects this
2. **Different parameter updates** - Algorithm prioritizes arr-h and heads
3. **Better visualization** - Plots show meaningful differences between groups

### After Convergence:
1. **Recession phi may be higher** - That's OK if it was unachievable anyway
2. **Arr-h and head phi should be lower** - Better match to achievable targets
3. **Overall calibration may improve** - If recession targets were unrealistic

## My Recommendation for local_run33

Based on your results:

**Start with Option 1 (Moderate Reduction):**

```python
# New weights in phi_factors.csv or parameterization:
3obg10:  0.05 → 0.015  # recession:pk4
22obg8:  0.05 → 0.015  # recession:pk4-conf-diff
23obg9:  0.05 → 0.010  # recession:pk4-spr-diff (worst performer, reduce more)

# Keep others the same:
arr-h groups: 0.05-0.30 (no change)
head groups: 0.30 (no change)
budget groups: 0.20 (no change)
```

**Why this approach:**
1. ✅ Recession contribution drops to ~25-30%
2. ✅ ARR-H and other groups become clearly visible
3. ✅ Still attempts to match recession (just with less influence)
4. ✅ Can tighten further if needed based on results
5. ✅ Easy to revert if recession performance degrades too much

## Monitoring After Weight Change

Run the phi analysis script after your next PEST++ run:

```bash
python scripts/f_analyze_phi.py local_run34
```

Look for:
- [ ] Recession contribution is 20-40% (good balance)
- [ ] ARR-H contribution is 2-5% (visible in plots)
- [ ] Recession phi reduction still >60% (not totally ignored)
- [ ] Arr-h and head phi reduction >99% (maintained excellence)
- [ ] Plots clearly show performance differences

## Summary Table

| Option | Recession Weight | Expected Recession Contrib | Visual Clarity | Calibration Balance |
|--------|-----------------|---------------------------|----------------|---------------------|
| Current | 0.05 | 96% | ❌ Poor | ⚠️ Recession-dominated |
| **Option 1** | **0.015** | **~30%** | **✅ Good** | **✅ Balanced** |
| Option 2 | 0.005 | ~10% | ✅ Excellent | ⚠️ May ignore recession |
| Option 3 | 0.005-0.05 | ~40% | ⚠️ Moderate | ✅ Selective |
| Option 4 | 0.001→0.02 | Staged | ✅ Excellent | ✅ Optimal (complex) |

**My vote: Option 1 with new weights of 0.010-0.015 for recession observations** 🎯
