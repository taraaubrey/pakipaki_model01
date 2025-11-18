# Phi Analysis - Consistent Color Scheme Applied to All Plots

## ✅ Complete Implementation

The grouped color scheme is now consistently applied across **ALL plots** in the phi analysis visualization.

## Color Scheme Reference

### 🔴 Red/Orange Shades - ARR Boundary Conditions
| Observation | Color | Hex Code | Description |
|-------------|-------|----------|-------------|
| arr-awq | Red | `#e74c3c` | Awanui stream array heads |
| arr-spq | Dark Red | `#c0392b` | Springs array heads |
| arr-confq | Orange | `#e67e22` | Confined area array heads |

### 🔵 Blue Shades - Head Observations
| Observation | Color | Hex Code | Description |
|-------------|-------|----------|-------------|
| arr-h | Blue | `#3498db` | Generic array heads |
| ts-heads | Dark Blue | `#2980b9` | Time series heads |

### 🟢 Green - Recession Observations
| Observation | Color | Hex Code | Description |
|-------------|-------|----------|-------------|
| recession | Green | `#2ecc71` | All recession rate observations |

### 🟣 Purple - Budget Observations
| Observation | Color | Hex Code | Description |
|-------------|-------|----------|-------------|
| budget | Purple | `#9b59b6` | All budget component observations |

## Plots Using Consistent Color Scheme

### ✅ Plot 1: Total Phi Evolution
- Black line with gray shading for ±1 std dev
- Red dashed lines for reinflation iterations

### ✅ Plot 2: ARR-H Observations Phi Evolution
**Colors:**
- arr-awq:awq → 🔴 Red
- arr-spq:spq → 🔴 Dark Red
- arr-confq:confq → 🟠 Orange
- arr-h → 🔵 Blue

### ✅ Plot 3: Head Observations Phi Evolution
**Colors:**
- ts-heads:pk4 → 🔵 Dark Blue
- ts-heads:pk4-conf-diff → 🔵 Dark Blue
- ts-heads:pk4-spr-diff → 🔵 Dark Blue

### ✅ Plot 4: Recession Observations Phi Evolution
**Colors:**
- recession:pk4 → 🟢 Green
- recession:pk4-conf-diff → 🟢 Green
- recession:pk4-spr-diff → 🟢 Green

### ✅ Plot 5: Normalized Phi Evolution (All Groups)
**Colors applied to all observation types:**
- ARR boundary conditions → Red/Orange shades
- Head observations → Blue shades
- Recession observations → Green
- Budget observations → Purple

### ✅ Plot 6: Final Phi Values (Non-Zero Weights Only)
**Horizontal bar chart with grouped colors**
- Shows only 12 active observations
- Same color scheme as other plots

### ✅ Plot 7: Phi Reduction Percentages (Non-Zero Weights Only)
**Vertical bar chart with grouped colors and legend**
- Shows only 12 active observations
- Dynamic legend showing only present observation types
- Same color scheme as other plots

## Visual Consistency Benefits

1. **Instant Recognition**: Same colors mean same observation types across all plots
2. **Easy Comparison**: Track specific observations across different visualizations
3. **Professional Appearance**: Consistent branding throughout the figure
4. **Better Storytelling**: Color groups tell the story of observation categories

## Example Interpretation

Looking at the plots for local_run33, you can now quickly identify:

- **Red/Orange lines/bars** = ARR boundary heads (performing well, ~99.6-99.97% reduction)
- **Blue lines/bars** = Head observations (performing excellently, ~99.98-100% reduction)
- **Green lines/bars** = Recession rates (struggling, ~72-78% reduction) ⚠️
- **Purple lines/bars** = Budget components (performing well, ~99.8-99.9% reduction)

This makes it immediately obvious that **green (recession) observations are the primary contributors to residual phi**, accounting for 96.48% of total residual.

## Usage

Same command, better visualizations:

```bash
# Default (uses MODEL_NAME from setup.py)
python scripts/f_analyze_phi.py

# For any specific run
python scripts/f_analyze_phi.py local_run34
```

## Complete Feature List

| Feature | Status |
|---------|--------|
| usecol appended to obs names | ✅ |
| Plot 2-4: Grouped color scheme | ✅ |
| Plot 5: Grouped color scheme | ✅ |
| Plot 6: Non-zero weights + colors | ✅ |
| Plot 7: Non-zero weights + colors | ✅ |
| Dynamic legend (only present colors) | ✅ |
| Reinflation markers (all plots) | ✅ |
| Consistent styling across plots | ✅ |

## Files

- **Script**: `scripts/f_analyze_phi.py` - Fully updated
- **Documentation**:
  - `scripts/README_f_analyze_phi.md` - Usage guide
  - `COLOR_SCHEME_COMPLETE.md` - This document
- **Example outputs**: `models/local_run33/local_run33_phi_analysis.png`

The phi analysis tool is now complete and ready for production use! 🎨📊
