# ARR-H Observation Group Analysis Summary
## PEST++ IES Run: local_run33

---

## Key Findings

### 1. **ARR-H Groups ARE Working, But Show Different Behavior Pattern**

The arr-h observation groups are actually achieving **very high phi reduction** (>99.6%), contrary to initial concerns. However, they show a distinct instability pattern compared to head observations.

### 2. **ARR-H Performance Summary**

| Group | Description | Weight | Initial Phi | Final Phi | Reduction % |
|-------|-------------|--------|-------------|-----------|-------------|
| **18obg4** | arr-awq (Awanui stream) | 0.05 | 4.97e+14 | 2.41e+11 | **99.95%** |
| **19obg5** | arr-spq (Springs) | 0.05 | 4.97e+14 | 1.52e+11 | **99.97%** |
| **greater_20obg6** | arr-confq (Confined) | 0.10 | 9.94e+14 | 3.63e+12 | **99.63%** |

**Comparison**: Head observations (1obg0, 2obg1, 13obg2) achieved 99.997% reduction, only marginally better.

---

## 3. **The Critical Issue: Iteration 4 Instability**

### Problem: Massive Phi Spike at Iteration 4

All three arr-h groups experienced catastrophic phi increases at **iteration 4**:

| Group | Iter 3 Phi | Iter 4 Phi | Increase Factor |
|-------|------------|------------|-----------------|
| **18obg4** (Awanui) | 1.10e+12 | 2.08e+17 | **×190,000** |
| **19obg5** (Spring) | 1.33e+13 | 8.69e+14 | **×65** |
| **greater_20obg6** (Confined) | 4.40e+11 | 1.89e+14 | **×430** |

This suggests that:
- The parameter ensemble at iteration 4 produced extremely poor matches to arr-h constraints
- Head observations were less sensitive to this parameter update
- **Iteration 4 may have been a re-inflation iteration** (check `local_run33.4.reinflate.pcs.csv`)

---

## 4. **Variability Across Ensemble**

The coefficient of variation (CV) reveals high uncertainty in arr-h observations:

### Coefficient of Variation (%) by Iteration

| Iteration | 18obg4 (Awanui) | 19obg5 (Spring) | greater_20obg6 (Confined) |
|-----------|-----------------|-----------------|---------------------------|
| 0 | 1008% | 1571% | 460% |
| 1 | 1136% | 842% | 367% |
| 2 | 863% | 1146% | 350% |
| 3 | 866% | 1536% | 369% |
| **4** | **1580%** | **1447%** | **674%** |
| 5 | 670% | 1563% | 652% |
| 6 | 1082% | 1039% | 651% |
| 7 | 692% | 1394% | 564% |
| 8 | 724% | 1031% | 628% |

**Interpretation**:
- Extremely high CV (>400%) indicates **large spread across realizations**
- Some realizations match arr-h observations well, others very poorly
- This suggests **parameter combinations have highly variable effects** on head arrays

---

## 5. **Comparison with Non-ARR-H Groups**

### Best vs Worst Performers (Final Reduction %)

**Top Performers:**
1. Head obs 13 (13obg2): 99.997% ⬅ **Non-ARR-H**
2. Head obs 1 (1obg0): 99.997% ⬅ **Non-ARR-H**
3. Head obs 2 (2obg1): 99.997% ⬅ **Non-ARR-H**
4. Head obs less_16 (less_16obg22): 99.984% ⬅ **Non-ARR-H**
5. **arr-spq (19obg5): 99.969%** ⬅ **ARR-H (BEST)**
6. **arr-awq (18obg4): 99.952%** ⬅ **ARR-H**

**Worst Performers:**
1. Recession rate 23 (23obg9): 71.94%
2. Recession rate 22 (22obg8): 78.30%
3. Recession rate 3 (3obg10): 78.30%
4. **arr-confq (greater_20obg6): 99.63%** ⬅ **ARR-H (WORST, but still 99.6%!)**

**Key Insight**: ARR-H groups rank in the middle-to-upper tier, NOT the worst performers. The recession rate observations (23obg9, 22obg8, 3obg10) are struggling far more.

---

## 6. **Contribution to Total Phi**

At final iteration (8), observation group contributions:

| Rank | Group | Contribution % | Type |
|------|-------|----------------|------|
| 1 | 23obg9 (Recession 23) | 37.98% | Non-ARR-H |
| 2 | 22obg8 (Recession 22) | 29.37% | Non-ARR-H |
| 3 | 3obg10 (Recession 3) | 29.13% | Non-ARR-H |
| 4 | less_6obg13 (Head) | 1.28% | Non-ARR-H |
| 5 | **greater_20obg6 (Confined arr-h)** | **0.98%** | **ARR-H** |
| 6 | greater_5obg12 (Head) | 0.94% | Non-ARR-H |
| 7 | less_16obg22 (Head) | 0.13% | Non-ARR-H |
| 8 | **18obg4 (Awanui arr-h)** | **0.07%** | **ARR-H** |
| 9 | **19obg5 (Spring arr-h)** | **0.04%** | **ARR-H** |

**Critical Finding**: ARR-H groups contribute **only 1.09% of total residual phi** at final iteration, meaning they are **well-matched** compared to recession rate observations.

---

## 7. **Why ARR-H Observations Appear Problematic (But Aren't)**

### Common Misconceptions:

❌ **"ARR-H isn't working"**
✅ **Reality**: ARR-H achieved 99.6-99.97% reduction, matching nearly as well as head obs

❌ **"ARR-H has the lowest phi"**
✅ **Reality**: ARR-H has the **smallest contribution to final phi** (1.09% total), meaning they're among the best-fitted

❌ **"ARR-H is unstable"**
✅ **Reality**: ARR-H shows instability at iteration 4, but **recovers**. This is likely due to re-inflation or aggressive parameter updates, not observation failure.

### Actual Issues:

1. **High ensemble variability (CV > 600%)**: Some parameter sets produce vastly different head arrays
2. **Iteration 4 spike**: Re-inflation or parameter update caused temporary poor fit
3. **Low weight (0.05-0.10)**: ARR-H observations have lower influence than head obs (0.20-0.30), so parameter updates prioritize head observations

---

## 8. **Recommendations**

### If You Want ARR-H to Have More Influence:

1. **Increase phi weights** in `phi_factors.csv`:
   - Current: 18obg4=0.05, 19obg5=0.05, greater_20obg6=0.10
   - Suggested: Increase to 0.10-0.20 if head array patterns are critical

2. **Investigate iteration 4**:
   - Check `local_run33.4.reinflate.pcs.csv` to see if re-inflation occurred
   - Review parameter updates between iterations 3→4
   - Consider adjusting lambda or re-inflation settings in PEST++ control file

3. **Examine observation weights**:
   - Current: Most arr-h obs have weight=0.0, only ~5% have weight=0.556
   - Consider increasing the proportion of non-zero weights for more spatial coverage

4. **Diagnose high ensemble spread**:
   - CV of 600-1500% suggests parameter combinations produce highly variable head fields
   - Check if this is physical (real parameter uncertainty) or numerical (model instability)
   - Consider tighter parameter bounds or prior constraints

### If Current Performance is Acceptable:

✅ **ARR-H observations are performing well** with 99.6-99.97% reduction
✅ They contribute minimal residual phi (1.09% total)
✅ The instability at iteration 4 was recovered by iteration 8

**Conclusion**: The arr-h observations are **not the problem**—they're actually among the best-performing groups. The recession rate observations (71-78% reduction) are the primary source of residual phi.

---

## 9. **Visual Analysis**

See generated plots:
- `models/local_run33/phi_analysis.png` - Overall phi evolution
- `models/local_run33/arr_h_detailed_analysis.png` - Detailed ARR-H comparison

Key visual findings:
- ARR-H and head obs follow similar downward trends (log scale)
- Iteration 4 spike visible in ARR-H groups
- Final phi values for ARR-H are orders of magnitude lower than recession obs
