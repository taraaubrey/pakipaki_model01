# Run 4 - PEST++ IES Diagnostic Analysis & Recommendations

**Date:** 2025-10-29
**Analysis:** PEST++ IES Record File (`models/local_run4/pest/master_ies/local_run4.rec`)
**Run Duration:** 12.7 hours (764.65 minutes)

---

## Executive Summary

Run 4 achieved **excellent overall performance** with a **92% reduction in phi** from initialization to best fit:
- **Initial Phi:** 3.24 × 10¹²
- **Best Phi:** 2.59 × 10¹¹ (iteration 12)
- **Final Phi:** 2.59 × 10¹¹

However, **prior-data conflict** exists with confined boundary observations, and **timeseries head observations dominate phi** (48%), indicating opportunities for improvement through parameter bound adjustments and phi factor rebalancing.

---

## Key Findings

### 1. Prior-Data Conflict (CRITICAL ISSUE)

**Status:** 2 observations in conflict with prior ensemble

| Observation | Target | Simulated Mean | Issue |
|------------|--------|----------------|-------|
| `confined_kper:2` | 0.0 m³/day | -9,250 m³/day | Boundary losing water |
| `confined_kper:4` | 0.0 m³/day | -7,835 m³/day | Boundary losing water |

**Interpretation:** The prior parameter ranges **cannot produce** the observed zero flux at the confined boundary. The model wants to drain water from the confined area, but observations indicate it should be neutral/balanced.

### 2. Observation Group Contributions (at Best Phi)

After phi factor weighting adjustment:

| Group | Observations | Phi Factor | % of Total Phi | Status |
|-------|--------------|------------|----------------|--------|
| **TS-HEADS** | 159 | 200 | **47.9%** | Too high |
| **RECESSION** | 3 | 50 | **26.2%** | Too high (only 3 obs!) |
| **GREATER_BUDGET** | 4 | 10 | **16.2%** | Conflicted |
| **LESS_BUDGET** | 8 | 10 | **9.0%** | Moderate |
| **HEAD-ARR** | 179,882 | 100 | **0.36%** | Well fit ✓ |
| **AWQ** | 121 | 1 | **0.37%** | Well fit ✓ |

### 3. Parameters at Bounds

Many critical parameters hitting bounds (potential for improvement):

| Parameter Group | Count | % at Bounds |
|----------------|-------|-------------|
| `ghb-conf-head-fngr` | 1,536 | **35%** |
| `ghb-aw-cond-fngr` | 430 | **33%** |
| `ghb-aw-head-fngr` | 430 | **28%** |

---

## Primary Issues

### Issue 1: Timeseries Heads Dominating Phi (48%)

**Problem:** TS-HEADS phi factor is 200 (highest) and contributes nearly half of total phi.

**Implications:**
- Parameter estimation overly focused on matching temporal dynamics
- May be pulling parameters away from fitting other observation types
- High phi factor relative to data content (159 observations)

**Root Causes:**
- Storage parameters (specific storage) may have insufficient range
- Recharge temporal patterns may not be flexible enough
- GHB head time series parameterization too rigid (temporal correlation constraints)

### Issue 2: Recession Observations Over-Weighted (26%)

**Problem:** Only 3 recession rate observations contributing 26% of phi (phi factor = 50)

**Implications:**
- Extreme leverage from very few observations
- Small errors in these 3 observations drive large parameter changes
- Parameter instability risk

### Issue 3: Confined Boundary Prior-Data Conflict

**Problem:** Structural/conceptual model issue - prior cannot produce observed fluxes

**Implications:**
- No combination of current parameters will eliminate this conflict
- IES will continue to struggle with this mismatch
- May bias other parameter estimates to compensate

**Potential Causes:**
- GHB conductance bounds too restrictive (can't achieve near-zero flux)
- Confined head parameterization insufficient (heads need to be higher)
- Boundary condition type may be inappropriate (consider no-flow or CHD instead)

---

## Recommendations

### Recommendation 1: Address Confined Boundary Prior-Data Conflict

**Priority:** CRITICAL
**Location:** `scripts/c_parameterization.py`

#### Option A: Increase GHB Conductance Bounds (PREFERRED)

- [ ] Increase `ghb-conf-cond` upper bound by 100-200%
- [ ] Review `ghb-conf-cond` geostatistical parameters (increase variance)
- [ ] Verify conductance parameterization allows sufficient spatial variability

```python
# In c_parameterization.py, look for ghb-conf-cond parameter setup
# Example modification:
# OLD: upper_bound = 100
# NEW: upper_bound = 300
```

#### Option B: Adjust Confined Head Parameterization

- [ ] Increase `ghb-conf-head` upper bound by 50%
- [ ] Increase `ghb-conf-head` spatial correlation length (allow smoother head field)
- [ ] Review temporal correlation - may be too restrictive

```python
# Allow higher heads to reduce hydraulic gradient
# Higher heads in confined boundary = less outflow
```

#### Option C: Reconsider Boundary Condition Type

- [ ] Investigate if confined boundary should be no-flow (if truly isolated)
- [ ] Consider constant head (CHD) if heads are well-known
- [ ] Review conceptual model with field data

**Verification:**
- [ ] Run prior ensemble only and check if flux observations can be matched
- [ ] Review `local_run5.0.pdc.csv` (next run) - should have 0 conflicts

---

### Recommendation 2: Rebalance Phi Factors

**Priority:** MODERATE
**Location:** `models/local_runX/pest/master_ies/phi_factors.csv`

#### Current Factors (Run 4):
```csv
head-arr,100
ts-heads,200
recession,50
AWq,1
budget,10
```

#### Recommended Factors (Run 5):

- [ ] Update phi_factors.csv with new values:

```csv
head-arr,100        # Keep - spatial constraints fitting well
ts-heads,100        # REDUCE from 200 - currently dominating phi
recession,25        # REDUCE from 50 - only 3 obs, too much leverage
AWq,1               # Keep - fitting well
budget,50           # INCREASE from 10 - need better flux match
```

**Rationale:**
- **TS-HEADS 200→100:** Reduce dominance, allow other obs types to influence fit
- **RECESSION 50→25:** Reduce leverage of 3 observations (or collect more recession data)
- **BUDGET 10→50:** Increase emphasis on resolving confined boundary flux issues

**Expected Impact:**
- More balanced phi contributions (no single group >35%)
- Improved budget fit (especially confined boundary)
- Slightly relaxed timeseries fit (acceptable trade-off)

---

### Recommendation 3: Widen Parameter Bounds

**Priority:** MODERATE
**Location:** `scripts/c_parameterization.py`

#### Parameters Hitting Bounds:

- [ ] **ghb-conf-head-fngr** (35% at bounds)
  - Increase upper bound by 0.5-1.0 m
  - Consider reducing lower bound if physically realistic

- [ ] **ghb-aw-head-fngr** (28% at bounds)
  - Review temporal head patterns - may need more flexibility
  - Increase variance in geostatistical parameters

- [ ] **ghb-aw-cond-fngr** (33% at bounds)
  - Increase upper bound by 50%
  - Check if conductance distribution is log-normal (should be)

- [ ] **Storage parameters** (stosslayer1-*)
  - Review specific storage bounds
  - Consider order-of-magnitude wider range

**Implementation:**
```python
# In c_parameterization.py, for each parameter:
# Check current bounds in pst.parameter_data
# Increase bounds by recommended percentages
# Ensure log-transformed parameters have appropriate bounds
```

---

### Recommendation 4: Investigate Temporal Dynamics

**Priority:** LOW-MODERATE
**Location:** `scripts/c_parameterization.py`, `scripts/e_pest_IES_HM.py`

Since TS-HEADS dominates phi, investigate:

- [ ] Review storage parameter correlation structure
  - Are specific storage values properly correlated with K?
  - Consider adding cross-correlation between npf-K and stos-S

- [ ] Examine recharge temporal correlation
  - May be too rigid (check temporal correlation length)
  - Consider allowing more inter-period variability

- [ ] Analyze GHB head time series parameterization
  - Review correlation structure across time steps
  - May need looser temporal correlation

- [ ] Identify which timesteps contribute most to TS-HEADS phi
  - Focus parameterization on poorly-fit periods
  - May reveal specific process deficiencies

---

### Recommendation 5: Add Regularization (Optional)

**Priority:** LOW
**Location:** `scripts/c_parameterization.py`

With 6,300+ parameters and 295 observations:

- [ ] Implement Tikhonov regularization
  - Add preferred value regularization toward geologically reasonable values
  - Weight regularization appropriately (don't overwhelm observations)

- [ ] Add parameter difference regularization
  - Prefer smooth parameter fields
  - Reduces oscillatory parameter patterns

**Example Implementation:**
```python
# In c_parameterization.py, after setting up prior:
# pst.add_regularization(...)
# Set preferred values based on field data or prior knowledge
# Adjust regularization weight (start with 0.01-0.1 × observation weight)
```

---

### Recommendation 6: Enhance Recession Observations (Data Collection)

**Priority:** LOW (Long-term)
**Location:** Field data collection

Current issue: Only 3 recession observations with 26% phi contribution

- [ ] Collect additional recession rate measurements
  - Multiple wells/springs across domain
  - Multiple recession events

- [ ] OR reduce recession phi factor to 10-15 (if data cannot be collected)

- [ ] Consider using recession as a regularization target rather than hard observation

---

## Implementation Checklist

### Immediate Actions (Run 5):

- [ ] Update `phi_factors.csv`:
  - TS-HEADS: 200 → 100
  - RECESSION: 50 → 25
  - BUDGET: 10 → 50

- [ ] Modify `c_parameterization.py` parameter bounds:
  - `ghb-conf-head`: +50% upper bound
  - `ghb-conf-cond`: +100% upper bound
  - `ghb-aw-cond`: +50% upper bound

- [ ] Verify changes by running prior ensemble only:
  ```bash
  # In d_pest_IES_prior.py, check prior-data conflict
  ```

### Short-term Actions (Run 6-7):

- [ ] Review storage parameter bounds (stosslayer1-*)
- [ ] Examine temporal correlation structures (recharge, GHB heads)
- [ ] Add Tikhonov regularization
- [ ] Verify truth data for confined boundary (is zero flux realistic?)
- [ ] Analyze per-timestep phi contributions for TS-HEADS

### Long-term Actions:

- [ ] Consider alternative boundary condition for confined area (no-flow or CHD)
- [ ] Collect additional recession rate observations
- [ ] Re-evaluate conceptual model if conflicts persist
- [ ] Investigate parameter correlation structure (K vs S)

---

## Expected Outcomes

With recommended changes (Run 5):

- **Phi Reduction:** Additional 30-50% reduction
  - Target: **1.0 - 1.5 × 10¹¹** (vs current 2.59 × 10¹¹)

- **Prior-Data Conflict:** Eliminated
  - Confined boundary observations should be fitable

- **Observation Group Balance:**
  - No single group >35% of phi
  - More distributed parameter sensitivity

- **Parameter Stability:**
  - <20% of parameters at bounds
  - Reduced oscillatory behavior

- **Overall Model Quality:**
  - Better representation of system dynamics
  - More reliable uncertainty quantification
  - Improved forecast confidence

---

## Analysis Notes

### Files Referenced:
- `models/local_run4/pest/master_ies/local_run4.rec` - Main record file
- `models/local_run4/pest/master_ies/local_run4.0.pdc.csv` - Prior-data conflict listing
- `models/local_run4/pest/master_ies/phi_factors.csv` - Phi factor weightings
- `truth/output.sample_heads.truth.csv` - Observation truth data

### Phi Progression:
```
Iteration:  0  1     2      3        4      5      6      ...  12     14
Phi:        3.24 2.68  1.01   1.02e+2  1.09   0.74   0.50  ...  0.26   1.34
            (×10¹²)                                               (best)
```

### Key Metrics:
- **Total Parameters:** 6,300+
- **Total Observations:** 295 (weighted)
- **Ensemble Size:** 200 realizations
- **IES Iterations:** 14 (noptmax = 11, continued with reinflation)
- **Acceptable Phi Factor:** 1.05

---

## Next Steps

1. Implement immediate actions checklist above
2. Run `c_parameterization.py` with updated bounds
3. Copy phi_factors.csv to new run directory (run5) with updated values
4. Execute `e_pest_IES_HM.py` for run5
5. Review `run5.rec` for:
   - Prior-data conflict resolution
   - Balanced observation group contributions
   - Improved phi reduction
6. Document results and iterate as needed

**Good luck with Run 5!**
