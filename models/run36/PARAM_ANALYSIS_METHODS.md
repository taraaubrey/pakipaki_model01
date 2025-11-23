# Parameter Analysis Methods

## Summary

This document describes how the mean and coefficient of variation (CV) are calculated in the combined parameter analysis script (`m_combined_param_analysis.py`).

## Calculations

### For Each Parameter Group (per iteration)

1. **Mean**: Arithmetic mean of all parameter values across all realizations
   ```
   mean = np.mean(group_values)
   ```
   where `group_values` contains all parameter values for that parameter group across all ensemble realizations.

2. **Standard Deviation**: Population standard deviation
   ```
   std = np.std(group_values)
   ```

3. **Coefficient of Variation (CV)**: Relative measure of variability expressed as a percentage
   ```
   CV = (std / mean) * 100
   ```
   CV is undefined (NaN) when mean = 0.

### Aggregation for Plots

For the evolution plots (B-E for mean, F-I for CV):

- **Mean by iteration**: For each subcategory, the mean values are averaged across all parameter groups within that subcategory
  ```python
  mean_by_iter = subcat_data.groupby('Iteration')['Mean'].mean()
  ```

- **CV by iteration**: For each subcategory, the CV values are averaged across all parameter groups within that subcategory
  ```python
  cv_by_iter = subcat_data.groupby('Iteration')['CV'].mean()
  ```

## Interpretation

### Mean Evolution (Plots B-E)
- Shows how the average parameter values change through IES iterations
- **Decrease**: Parameters are being reduced during history matching
- **Increase**: Parameters are being increased during history matching
- **Stability**: Parameters have converged to stable values
- **Jumps at reinflation**: Parameter ensemble is reinflated to increase variability

### CV Evolution (Plots F-I)
- Shows how the parameter uncertainty (spread) changes through iterations
- **CV decrease**: Uncertainty reduction - the ensemble is collapsing toward preferred values
- **CV increase**: Uncertainty increase - typically occurs at reinflation iterations
- **High CV (>100%)**: Large uncertainty in parameter values
- **Low CV (<10%)**: Parameters are well constrained

### Units
- **Mean**: Same units as the parameter (e.g., m/d for hydraulic conductivity)
- **CV**: Percentage (%)

## Notes

- Log-scale is used for Hydraulic Conductivity and GHB Conductance mean plots due to their wide value ranges
- Reinflation iterations are marked with red dashed vertical lines
- Each subcategory within a parameter type uses a distinct color from a gradient
