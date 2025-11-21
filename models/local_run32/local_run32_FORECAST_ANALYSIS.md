# Forecast Analysis Summary: local_run32

## Overview

Total forecasts: 44

Iterations analyzed: 17

Last iteration: 16

Reinflation iterations: 2, 4, 6, 8, 9, 10, 12, 15, 16, 18, 20, 22, 24, 26, 28, 30, 32, 35, 37, 38, 40, 42, 43, 45, 48, 50, 53, 54, 55, 59, 62, 64, 66, 68, 70, 73, 74, 75, 77, 79, 81, 82, 84, 86, 88, 89, 93, 96, 97, 101, 103, 105, 106, 108, 111, 112, 114, 116, 117, 118, 120, 121, 123, 124, 125, 128, 132, 133, 135, 136, 138, 139, 142, 144, 146, 148, 150, 151, 153, 155, 156, 157, 159, 160, 163, 164, 165, 167, 172, 173, 175, 178, 179, 181, 182, 185, 186, 188, 189, 190, 191, 194, 195, 196, 198, 201, 203, 207, 208, 209, 213, 215, 218, 220, 222, 224, 227, 230, 233, 234, 237, 239, 240, 243, 246, 247

## Forecast Categories

| category       |   Count |   Avg Std Reduction % |   Avg CV Reduction % |   Truth in Post Range |
|:---------------|--------:|----------------------:|---------------------:|----------------------:|
| Budget-GHB     |       8 |                 81.27 |                17.53 |                     4 |
| Budget-Inflow  |       4 |                nan    |               nan    |                     4 |
| Budget-SW      |      12 |                 92.39 |               -11.87 |                    12 |
| Budget-Storage |       4 |                 60.84 |                -4.09 |                     4 |

## Top Uncertainty Reductions

| component      | period   |        prior_std |         post_std |   std_reduction_pct | truth_in_post   |
|:---------------|:---------|-----------------:|-----------------:|--------------------:|:----------------|
| Spring GHB     | kper1    |      1.5125e+06  |  20589           |             98.6387 | True            |
| Spring GHB     | kper2    | 514631           |  17039.2         |             96.6891 | True            |
| Awanui Stream  | kper1    |      3.53327e+07 |      1.43311e+06 |             95.9439 | True            |
| Confined GHB   | kper1    |      3.54405e+07 |      1.44722e+06 |             95.9165 | False           |
| Total SW       | kper1    |      3.54405e+07 |      1.44722e+06 |             95.9165 | True            |
| Confined GHB   | kper2    |      2.19615e+07 | 951895           |             95.6656 | False           |
| Awanui Stream  | kper2    |      2.19166e+07 | 958316           |             95.6274 | True            |
| Total SW       | kper2    |      2.19988e+07 | 968063           |             95.5995 | True            |
| Poukawa Stream | kper2    |      1.30296e+06 | 113025           |             91.3255 | True            |
| Poukawa Stream | kper1    |      1.93653e+06 | 170346           |             91.2036 | True            |

## Largest Forecast Mean Changes

| component     | period   |   prior_mean |         post_mean |   mean_change |   mean_change_pct |
|:--------------|:---------|-------------:|------------------:|--------------:|------------------:|
| Total SW      | kper3    | -1.73026e+07 |      -1.19973e+06 |   1.61029e+07 |          -93.0662 |
| Confined GHB  | kper3    |  1.73007e+07 |       1.19778e+06 |  -1.61029e+07 |          -93.0767 |
| Awanui Stream | kper3    | -1.72948e+07 |      -1.19845e+06 |   1.60964e+07 |          -93.0705 |
| Total SW      | kper4    | -1.40504e+07 | -996825           |   1.30536e+07 |          -92.9054 |
| Awanui Stream | kper4    | -1.40443e+07 | -995443           |   1.30488e+07 |          -92.9121 |
| Confined GHB  | kper4    |  1.40162e+07 |  981308           |  -1.30349e+07 |          -92.9988 |
| Confined GHB  | kper1    |  9.20609e+06 |  493788           |  -8.7123e+06  |          -94.6363 |
| Total SW      | kper1    | -9.20804e+06 | -495745           |   8.7123e+06  |          -94.6162 |
| Awanui Stream | kper1    | -8.69401e+06 | -467759           |   8.22625e+06 |          -94.6197 |
| Total SW      | kper2    | -5.80244e+06 | -347666           |   5.45477e+06 |          -94.0083 |

## Constraint Forecasts

| component    | period   | constraint_type   |   truth |        post_mean |         post_std |
|:-------------|:---------|:------------------|--------:|-----------------:|-----------------:|
| Confined GHB | kper1    | greater_than      |       0 | 493788           |      1.44722e+06 |
| Confined GHB | kper2    | greater_than      |       0 | 319875           | 951895           |
| Confined GHB | kper3    | greater_than      |       0 |      1.19778e+06 |      7.1607e+06  |
| Confined GHB | kper4    | greater_than      |       0 | 981308           |      5.45111e+06 |

