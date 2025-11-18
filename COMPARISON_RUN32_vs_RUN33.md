# Comparison: local_run32 vs local_run33

## Overview

| Metric | local_run32 | local_run33 | Change |
|--------|-------------|-------------|--------|
| **Total Iterations** | 16 | 8 | 🔽 Half as many |
| **Initial Phi (mean)** | 4.09e+11 | 1.94e+16 | ⬆️ 47,500× higher |
| **Final Phi (mean)** | 1.82e+09 | 3.70e+14 | ⬆️ 203,700× higher |
| **Overall Reduction** | 99.56% | 98.10% | 🔽 Slightly worse |
| **Reinflation Iterations** | 4, 10 | 4 | Same pattern |

## Key Differences

### 1. Scale of Problem

**local_run32**: Started with moderate initial phi (4.09e+11)
- More "reasonable" starting point
- Converged to very low final phi (1.82e+09)

**local_run33**: Started with very high initial phi (1.94e+16)
- 47,500× larger initial mismatch
- Still has higher residual despite good reduction
- Suggests different parameterization or truth data

### 2. Observation Weight Differences

| Observation | run32 Weight | run33 Weight | Impact |
|-------------|--------------|--------------|--------|
| **arr-h (less_16obg22)** | **0.40** | **0.30** | Run32 prioritized arr-h more |
| arr-awq, arr-spq | 0.05 | 0.05 | Same |
| arr-confq | 0.10 | 0.10 | Same |
| ts-heads | 0.30 | 0.30 | Same |
| recession | 0.05 | 0.05 | Same |
| budget | 0.20 | 0.20 | Same |

**Key finding**: Run32 had higher weight (0.40 vs 0.30) on arr-h observations!

## Performance by Observation Type

### ARR-H Observations

| Observation | run32 Reduction | run33 Reduction | run32 Final Contrib | run33 Final Contrib |
|-------------|-----------------|-----------------|---------------------|---------------------|
| **arr-h** | **99.9996%** ✅ | **99.98%** ✅ | 0.02% | 0.13% |
| arr-awq:awq | 99.52% | 99.95% ✅ | 2.66% | 0.07% |
| arr-spq:spq | 99.89% | 99.97% ✅ | 0.59% | 0.04% |
| arr-confq:confq | 99.97% | 99.63% | 0.29% | 0.98% |
| **Total ARR-H** | **3.56%** | **1.21%** | **Better in run33** ✅ |

**Analysis:**
- ✅ **run33 has better ARR-H performance** (1.21% vs 3.56% contribution)
- arr-h (less_16obg22) performed excellently in both runs
- run33 improved arr-awq and arr-spq significantly

### Head Observations

| Observation | run32 Reduction | run33 Reduction | run32 Final Contrib | run33 Final Contrib |
|-------------|-----------------|-----------------|---------------------|---------------------|
| ts-heads:pk4 | 99.998% | 99.997% | 0.07% | 0.03% |
| ts-heads:pk4-conf-diff | 99.998% | 99.997% | 0.07% | 0.03% |
| ts-heads:pk4-spr-diff | 99.998% | 99.997% | 0.07% | 0.03% |
| **Total Head** | **0.20%** | **0.08%** | **Better in run33** ✅ |

**Analysis:**
- ✅ **run33 has slightly better head performance**
- Both runs achieve >99.99% reduction (excellent)

### Recession Observations

| Observation | run32 Reduction | run33 Reduction | run32 Final Contrib | run33 Final Contrib |
|-------------|-----------------|-----------------|---------------------|---------------------|
| recession:pk4 | **95.27%** | **78.30%** | 25.87% | 29.13% |
| recession:pk4-conf-diff | **95.27%** | **78.30%** | 25.98% | 29.37% |
| recession:pk4-spr-diff | **92.76%** | **71.94%** | 39.66% | 37.98% |
| **Total Recession** | **91.44%** | **96.48%** | **Much worse in run33** ⚠️ |

**Analysis:**
- ⚠️ **run33 has significantly worse recession performance**
- run32: 92-95% reduction (acceptable)
- run33: 72-78% reduction (poor)
- run33 recession dominates residual phi even more (96% vs 91%)

### Budget Observations

| Observation | run32 Reduction | run33 Reduction | run32 Final Contrib | run33 Final Contrib |
|-------------|-----------------|-----------------|---------------------|---------------------|
| budget:confined | 99.96% | 99.83% | 0.95% | 0.94% |
| budget:inflow (flux) | 99.82% | 99.76% | 3.86% | 1.28% |
| **Total Budget/Flux** | **4.81%** | **2.22%** | **Better in run33** ✅ |

**Analysis:**
- ✅ **run33 has better budget performance**

## Final Phi Contribution Comparison

### run32 (Total residual: 1.82e+09)
```
Recession:  91.44% 🟢 (acceptable, 92-95% reduction)
Flux:        3.86% 🟣
ARR-H:       3.56% 🔴
Budget:      0.95% 🟣
Head:        0.20% 🔵
```

### run33 (Total residual: 3.70e+14)
```
Recession:  96.48% 🟢 (poor, 72-78% reduction) ⚠️
ARR-H:       1.21% 🔴
Flux:        1.28% 🟣
Budget:      0.94% 🟣
Head:        0.08% 🔵
```

## Key Insights

### 1. What Improved in run33 ✅
- ARR-H observations: 3.56% → 1.21% (better match)
- Head observations: 0.20% → 0.08% (better match)
- Budget observations: Similar performance
- Fewer iterations needed (16 → 8)

### 2. What Got Worse in run33 ⚠️
- **Recession observations**: 91.44% → 96.48% of total phi
- **Recession reduction**: 92-95% → 72-78%
- Much higher absolute phi values (scaling issue?)
- Recession now completely dominates the problem

### 3. Why the Difference?

**Possible explanations:**

1. **Different initial conditions or truth data**
   - Initial phi 47,500× higher suggests very different starting point
   - May have changed observation values or truth targets

2. **Recession observations became harder to match**
   - Reduction dropped from 92-95% to 72-78%
   - Something changed in recession calculation or targets
   - May need to review recession rate observation weights

3. **ARR-H improvement trade-off**
   - Better ARR-H performance might have come at cost of recession
   - Parameters may have shifted to favor spatial heads over temporal dynamics

## Recommendations Based on Comparison

### For run33:

1. **Reduce recession weights** (as discussed earlier)
   ```
   Current: 0.05
   Suggested: 0.015
   Reason: Recession dominates even more in run33 (96% vs 91% in run32)
   ```

2. **Investigate why recession got worse**
   - Check if recession rate calculation changed
   - Verify recession truth data is reasonable
   - Compare parameter values between run32 and run33
   - Check storage properties (specific storage, specific yield)

3. **Consider run32's arr-h weight**
   - run32 used 0.40 for arr-h (vs 0.30 in run33)
   - This gave excellent arr-h performance (99.9996% reduction)
   - May want to try run33 with 0.40 arr-h weight

### For future runs:

1. **Use run32 as baseline for comparison**
   - It achieved better balance (91% recession vs 96%)
   - Better recession reduction (92-95% vs 72-78%)
   - Lower overall phi (though on different scale)

2. **Monitor recession-to-ARR-H trade-off**
   - Improving one shouldn't degrade the other this much
   - May need to review parameter constraints

3. **Standardize initial phi scale**
   - Investigate why run33 started 47,500× higher
   - Ensure consistent normalization between runs

## Summary Table

| Aspect | Winner | Notes |
|--------|--------|-------|
| ARR-H Performance | **run33** ✅ | 1.21% vs 3.56% contribution |
| Head Performance | **run33** ✅ | 0.08% vs 0.20% contribution |
| Recession Performance | **run32** ✅ | 91.44% vs 96.48% contribution; 92-95% vs 72-78% reduction |
| Budget Performance | **run33** ✅ | Slightly better |
| Overall Balance | **run32** ✅ | Better distribution across observation types |
| Efficiency | **run33** ✅ | 8 iterations vs 16 |
| **Overall Recommendation** | **run32** | Better recession performance and balance |

## Action Items

- [ ] Review what changed between run32 and run33 parameterization
- [ ] Investigate recession observation calculation/targets
- [ ] Consider using run32's weights (especially 0.40 for arr-h)
- [ ] Reduce recession weights in run33 to 0.015 for better visualization
- [ ] Compare parameter values (K, S, recharge) between converged runs
