# Final Phi Analysis Script Updates

## Changes Implemented

### 1. ✅ Final Phi Plot - Now Shows Only Non-Zero Weighted Observations

**Plot 6: "Final Phi Values"** has been updated:
- Filters out all observations with weight=0
- Title updated to: "Final Phi Values (Non-Zero Weights Only)"
- Shows only the 12 active observations instead of 25
- Makes the plot cleaner and more relevant to calibration

### 2. ✅ New Color Scheme - Grouped by Observation Subcategory

**Previous scheme:** Generic colors by type (ARR-H=red, Head=blue, Recession=green)

**New scheme:** Specific colors by observation subcategory

#### Color Assignments:

| Observation Type | Color | Hex Code | Usage |
|-----------------|-------|----------|-------|
| **arr-awq** | Red | `#e74c3c` | Awanui stream array heads |
| **arr-spq** | Dark Red | `#c0392b` | Springs array heads |
| **arr-confq** | Orange | `#e67e22` | Confined area array heads |
| **arr-h** | Blue | `#3498db` | Generic array heads |
| **ts-heads** | Dark Blue | `#2980b9` | Time series heads |
| **recession** | Green | `#2ecc71` | Recession rate observations |
| **budget** | Purple | `#9b59b6` | Budget observations |

#### Grouping Logic:

1. **ARR stream/spring group** (arr-awq, arr-spq, arr-confq):
   - Shades of red/orange
   - Represents boundary condition heads
   - Differentiated by specific boundary type

2. **ARR-H and ts-heads group** (arr-h, ts-heads):
   - Shades of blue
   - Represents head observations
   - Differentiated by whether they're arrays or time series

3. **Recession group**:
   - Green
   - All recession observations

4. **Budget group**:
   - Purple
   - Budget component observations

### Affected Plots

The new color scheme is applied to:
- **Plot 6**: Final Phi Values (horizontal bar chart)
- **Plot 7**: Phi Reduction Percentages (vertical bar chart)

Both plots now:
- Show only non-zero weighted observations
- Use the grouped color scheme
- Include a legend showing only the observation types present in the data

### Example for local_run33

The updated plots now show **12 observations** (instead of 25), with colors:

| Observation | Color | Interpretation |
|-------------|-------|----------------|
| ts-heads:pk4 | Dark Blue | Head time series at pk4 |
| ts-heads:pk4-conf-diff | Dark Blue | Head diff from confined |
| ts-heads:pk4-spr-diff | Dark Blue | Head diff from spring |
| arr-h | Blue | Array head observations |
| arr-spq:spq | Dark Red | Springs array heads |
| arr-awq:awq | Red | Awanui stream array heads |
| arr-confq:confq | Orange | Confined area array heads |
| budget:confined | Purple | Budget for confined area |
| budget:inflow | Purple | Budget inflow component |
| recession:pk4 | Green | Recession at pk4 |
| recession:pk4-conf-diff | Green | Recession for confined diff |
| recession:pk4-spr-diff | Green | Recession for spring diff |

### Visual Benefits

1. **Easier to identify observation subcategories** at a glance
2. **Related observations are visually grouped** by color
3. **Cleaner plots** with only active (non-zero weighted) observations
4. **Better for comparing** similar observation types

### Usage

Same as before - no changes to command line interface:

```bash
# Default (uses MODEL_NAME from setup.py)
python scripts/f_analyze_phi.py

# For any specific run
python scripts/f_analyze_phi.py local_run34
```

### Complete Update Summary

| Feature | Status | Description |
|---------|--------|-------------|
| usecol appended to names | ✅ | Shows observation specifics (e.g., ts-heads:pk4) |
| Plot 7: Non-zero weights only | ✅ | Bar graph filters weight=0 obs |
| Plot 6: Non-zero weights only | ✅ | Final phi plot filters weight=0 obs |
| Grouped color scheme | ✅ | Colors group related observation types |
| Dynamic legend | ✅ | Shows only colors present in data |

### Files Updated

- `scripts/f_analyze_phi.py` - Main script with all updates
- `models/local_run33/local_run33_phi_analysis.png` - Regenerated with new colors
- `models/local_run33/local_run33_phi_summary.csv` - Updated with usecol names
- `models/local_run33/local_run33_PHI_ANALYSIS.md` - Updated report

### Next Steps

The script is now fully updated and ready to use for:
- Current run: local_run33
- Future runs: local_run34, local_run35, etc.

Just run the script after completing `e_pest_IES_HM.py` to get comprehensive phi analysis with the improved visualization scheme!

## Color Scheme Reference

Quick reference for interpreting the plots:

- 🔴 **Red shades** = ARR boundary conditions (awq, spq, confq)
- 🔵 **Blue shades** = Head observations (arr-h, ts-heads)
- 🟢 **Green** = Recession rates
- 🟣 **Purple** = Budget components
