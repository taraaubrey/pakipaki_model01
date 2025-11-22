# Parameter Analysis Summary: local_run34

## Overview

Total iterations completed: 20

Total parameter groups analyzed: 28

Parameter categories: GHB Head, GHB Conductance, Hydraulic Conductivity, Storage, Recharge

Initial phi (mean): 4.04e+15

Final phi (mean): 3.43e+14

Phi reduction: 91.50%

## Parameter Evolution by Category

### GHB Head

| Subcategory       | Pname             |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:------------------|:------------------|--------------:|-----------------:|---------------------:|---------------------:|
| GHB-Confined-head | ghbconf-head-hg   |   8351.77     |        98.542    |            0.0820801 |             0        |
| GHB-Spring-head   | ghbspring-head-hg |   5029.64     |        96.6573   |            0.5       |             0        |
| GHB-Awanui-head   | ghbaw-head-hg     |   1446.83     |        81.7982   |            0.0563246 |             0.199523 |
| GHB-Confined-head | ghbconf-head-cn   |    374.545    |        78.9794   |            0.4       |             0        |
| GHB-Spring-head   | ghbspring-head-cn |    347.59     |        77.6581   |            0         |             0        |
| GHB-Awanui-head   | ghbaw-head-cn     |    243.815    |        70.9479   |            0         |             0.4      |
| GHB-Poukawa-head  | ghbpw-head-hg     |   2076.73     |        36.9901   |            0.208     |             0        |
| GHB-Poukawa-head  | ghbpw-head-cn     |      0.670818 |         0.666348 |            0         |             0        |

### GHB Conductance

| Subcategory       | Pname             |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:------------------|:------------------|--------------:|-----------------:|---------------------:|---------------------:|
| GHB-Awanui-cond   | ghbaw-cond-cn     |     8407      |          98.7967 |             0        |             0.4      |
| GHB-Confined-cond | ghbconf-cond-gr   |     5583.97   |          98.6213 |             0.480956 |             0        |
| GHB-Confined-cond | ghbconf-cond-cn   |      746.26   |          88.2528 |             0.8      |             0        |
| GHB-Awanui-cond   | ghbaw-cond-gr     |     1954.31   |          82.99   |             0.133652 |             0.145107 |
| GHB-Poukawa-cond  | ghbpw-cond-gr     |     1019.38   |          77.2771 |             0.336    |             0.032    |
| GHB-Poukawa-cond  | ghbpw-cond-cn     |      167.618  |         -47.8886 |             0        |             0        |
| GHB-Spring-cond   | ghbspring-cond-gr |      222.252  |         -66.3517 |             0.266667 |             0.3      |
| GHB-Spring-cond   | ghbspring-cond-cn |       77.8793 |        -352.065  |             0        |             0        |

### Hydraulic Conductivity

| Subcategory   | Pname         |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:--------------|:--------------|--------------:|-----------------:|---------------------:|---------------------:|
| K-constant    | npfklayer1-cn |      3020.55  |          96.8426 |             0        |             3.2      |
| K-pilot       | npfklayer1-pp |      2583.69  |          95.7878 |             0        |             1.07586  |
| K-grid        | npfklayer1-gr |       805.601 |          30.882  |             0.198953 |             0.102618 |

### Storage

| Subcategory   | Pname          |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:--------------|:---------------|--------------:|-----------------:|---------------------:|---------------------:|
| S-constant    | stosslayer1-cn |       3744.05 |          97.268  |            0.4       |            0         |
| S-pilot       | stosslayer1-pp |       1761.2  |          84.8686 |            0.510345  |            0.0551724 |
| S-grid        | stosslayer1-gr |       1305.08 |          64.4561 |            0.0837696 |            0.131937  |

### Recharge

| Subcategory   | Pname    |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:--------------|:---------|--------------:|-----------------:|---------------------:|---------------------:|
| RCH-steady    | rchss-cn |      2234.68  |          95.7222 |            0.8       |            0         |
| RCH-transient | rchtr-pp |      1065.48  |          90.2883 |            0.717241  |            0         |
| RCH-steady    | rchss-gr |      1920.53  |          89.4995 |            0.270157  |            0.0712042 |
| RCH-steady    | rchss-pp |      1645.51  |          85.8937 |            0.317241  |            0         |
| RCH-transient | rchtr-gr |      1139.57  |          78.7294 |            0.0460733 |            0.272251  |
| RCH-transient | rchtr-cn |       130.178 |        -231.371  |            0         |            0         |

## Boundary Hitting Diagnostics

Parameters at bounds may indicate over-parameterization or improper bounds.

## Uncertainty Reduction Summary

| Category               |     mean |         min |     max |
|:-----------------------|---------:|------------:|--------:|
| GHB Conductance        | -2.54597 | -352.065    | 98.7967 |
| GHB Head               | 67.7799  |    0.666348 | 98.542  |
| Hydraulic Conductivity | 74.5041  |   30.882    | 96.8426 |
| Recharge               | 34.7937  | -231.371    | 95.7222 |
| Storage                | 82.1976  |   64.4561   | 97.268  |

