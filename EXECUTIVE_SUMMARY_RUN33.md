# Executive Summary: local_run33 PEST++ IES Analysis

**Date**: 2025-11-18
**Iterations**: 8 (with reinflation at iteration 4)
**Realizations**: 250 initial → 219 final

---

## 🎯 Key Findings

### 1. Overall Performance: Good But Unbalanced

| Metric | Value |
|--------|-------|
| **Overall Phi Reduction** | 98.10% |
| **Initial Phi (mean)** | 1.94e+16 |
| **Final Phi (mean)** | 3.70e+14 |
| **Iterations** | 8 |

**Verdict**: ✅ The model converges well, but one observation type dominates the residual.

---

## 📊 Observation Group Performance

### Excellent Performers (>99% reduction) ✅

| Observation Type | Reduction | Final Contribution |
|------------------|-----------|---------------------|
| **ts-heads:pk4** | 99.997% | 0.03% |
| **ts-heads:pk4-conf-diff** | 99.997% | 0.03% |
| **ts-heads:pk4-spr-diff** | 99.997% | 0.03% |
| **arr-h** | 99.98% | 0.13% |
| **arr-spq:spq** | 99.97% | 0.04% |
| **arr-awq:awq** | 99.95% | 0.07% |
| **budget:inflow** | 99.76% | 1.28% |
| **budget:confined** | 99.83% | 0.94% |
| **arr-confq:confq** | 99.63% | 0.98% |

**Total contribution from excellent performers**: 3.52%

### Poor Performers (<80% reduction) ⚠️

| Observation Type | Reduction | Final Contribution |
|------------------|-----------|---------------------|
| **recession:pk4-spr-diff** | 71.94% | 37.98% |
| **recession:pk4** | 78.30% | 29.13% |
| **recession:pk4-conf-diff** | 78.30% | 29.37% |

**Total contribution from recession**: 96.48% 🔴

---

## ⚠️ Problem #1: Recession Observations Dominate

**Issue**: Recession observations account for 96.48% of total residual phi, making other observation groups invisible in visualizations.

**Impact**:
- Hard to see performance of arr-h, heads, and budget observations
- Gives false impression that arr-h is performing poorly (it's actually excellent!)
- Calibration may be over-focusing on difficult-to-match recession rates

**Recommended Solution**:
```python
# Reduce recession weights from 0.05 to 0.015
# Expected result:
# - Recession contribution: ~30% (down from 96%)
# - ARR-H contribution: ~4% (up from 1.2%)
# - Better visualization without ignoring recession
```

**See**: `WEIGHT_ADJUSTMENT_RECOMMENDATIONS.md` for detailed options

---

## ⚠️ Problem #2: 39% of Realizations Are Stuck

**Issue**: 86 out of 219 realizations (39.3%) show minimal improvement:
- Less than 50% phi reduction
- Some actually got WORSE (negative reduction)
- Worst case: Realization 226 got **1199% worse** (phi increased 13×)

**Root Cause**: Low ensemble variance
- Stuck realizations receive **3× smaller parameter updates** than good ones
  - Stuck: 20,029 mean absolute change
  - Good: 58,950 mean absolute change
  - Ratio: 0.34×
- Parameters are NOT hitting bounds (only 1.3% at bounds)
- Algorithm can't find good directions for these realizations

**Recommended Solutions**:

1. **Increase ensemble diversity** (Priority #1):
   ```python
   # In d_pest_IES_prior.py
   prior_std_multiplier = 2.0  # Double the prior variance
   ```

2. **Increase ensemble size**:
   ```python
   NREALS_PRIOR = 500  # Up from 250
   ```

3. **Use subset selection**:
   ```python
   # In PEST++ control:
   ies_subset_size(100)  # Use best 100 realizations only
   ```

4. **More conservative lambda**:
   ```python
   ies_lambda_mults(0.5, 1.0, 2.0)  # Less aggressive
   ```

**See**: `STUCK_REALIZATIONS_DIAGNOSIS.md` for full analysis

---

## 📈 Comparison with local_run32

| Aspect | run32 | run33 | Winner |
|--------|-------|-------|--------|
| **Recession Performance** | 91.44% contrib | 96.48% contrib | run32 ✅ |
| **Recession Reduction** | 92-95% | 72-78% | run32 ✅ |
| **ARR-H Performance** | 3.56% contrib | 1.21% contrib | run33 ✅ |
| **Head Performance** | 0.20% contrib | 0.08% contrib | run33 ✅ |
| **Overall Balance** | Better | Recession-dominated | run32 ✅ |
| **Iterations** | 16 | 8 | run33 ✅ |

**Key Insight**: run32 had better recession performance and more balanced calibration. The change between runs seems to have improved spatial head matching at the expense of temporal dynamics (recession).

**See**: `COMPARISON_RUN32_vs_RUN33.md` for detailed comparison

---

## 🎯 Action Items for Next Run (local_run34)

### Priority 1: Fix Stuck Realizations
- [ ] Increase prior parameter variance by 2×
- [ ] Consider increasing ensemble size to 500

### Priority 2: Improve Visualization
- [ ] Reduce recession weights from 0.05 to 0.015

### Priority 3: Investigate Recession Issue
- [ ] Review what changed between run32 and run33
- [ ] Check recession rate calculation/targets
- [ ] Compare parameter values (K, S, recharge) between runs
- [ ] Consider run32's observation weights (arr-h was 0.40 vs 0.30)

### Priority 4: Refine PEST++ Settings
- [ ] Try more conservative lambda schedule: `ies_lambda_mults(0.5, 1.0, 2.0)`
- [ ] Consider subset selection: `ies_subset_size(100)`

---

## 📁 Generated Files

### Analysis Outputs:
- `models/local_run33/local_run33_phi_analysis.png` - 7-panel visualization with consistent color scheme
- `models/local_run33/local_run33_phi_summary.csv` - Performance metrics by observation type
- `models/local_run33/local_run33_PHI_ANALYSIS.md` - Detailed phi analysis report
- `models/local_run33/stuck_realizations.csv` - List of 86 stuck realizations
- `models/local_run33/realization_performance.csv` - All 219 realizations with phi reduction

### Documentation:
- `STUCK_REALIZATIONS_DIAGNOSIS.md` - Root cause analysis and solutions for stuck realizations
- `WEIGHT_ADJUSTMENT_RECOMMENDATIONS.md` - Detailed guide for adjusting observation weights
- `COMPARISON_RUN32_vs_RUN33.md` - Side-by-side comparison of the two runs
- `COLOR_SCHEME_COMPLETE.md` - Color scheme documentation for plots
- `EXECUTIVE_SUMMARY_RUN33.md` - This document

### Tools:
- `scripts/f_analyze_phi.py` - Phi analysis tool (works for any run)
- `analyze_stuck_realizations.py` - Diagnostic tool for parameter movement analysis

---

## 🚀 Quick Commands

```bash
# Run phi analysis on any model
python scripts/f_analyze_phi.py local_run34

# Analyze stuck realizations (update model_name in script)
python analyze_stuck_realizations.py

# Start a new run with recommended changes
# 1. Edit scripts/setup.py: MODEL_NAME = 'local_run34'
# 2. Edit scripts/d_pest_IES_prior.py: increase prior variance
# 3. Edit phi_factors.csv: reduce recession weights to 0.015
python scripts/a_model_setup.py
python scripts/c_parameterization.py
python scripts/d_pest_IES_prior.py
python scripts/e_pest_IES_HM.py
python scripts/f_analyze_phi.py
```

---

## 📝 Summary

**What's Working**:
- arr-h, ts-heads, and budget observations perform excellently (>99% reduction)
- Overall convergence is good (98% phi reduction)
- Fewer iterations than run32 (8 vs 16)

**What Needs Attention**:
- Recession observations dominate visualization (96% of residual)
- 39% of realizations stuck due to low ensemble variance
- Recession matching worse than run32 (72-78% vs 92-95% reduction)

**Recommended Focus for Next Run**:
1. Increase ensemble diversity (2× prior variance, 500 realizations)
2. Reduce recession weights to 0.015 for better balance
3. Investigate why recession performance degraded from run32
4. Consider using subset selection to improve efficiency

With these changes, you should see:
- <10% stuck realizations (down from 39%)
- Better visualization of all observation types
- Maintained excellent performance on arr-h and heads
- Potential improvement in recession matching

---

**Status**: ✅ Analysis complete and ready for next model run
