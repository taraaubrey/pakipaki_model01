# Parameter Analysis Summary: local_run32

## Overview

Total iterations completed: 16

Total parameter groups analyzed: 22

Parameter categories: GHB Head, GHB Conductance, Hydraulic Conductivity, Storage

Initial phi (mean): 4.09e+11

Final phi (mean): 1.82e+09

Phi reduction: 99.56%

## Parameter Evolution by Category

### GHB Head

| Subcategory       | Pname             |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:------------------|:------------------|--------------:|-----------------:|---------------------:|---------------------:|
| GHB-Spring-head   | ghbspring-head-cn |     30499.8   |          99.6891 |            0         |            0.97561   |
| GHB-Confined-head | ghbconf-head-cn   |     24487.8   |          99.5841 |            0         |            0         |
| GHB-Poukawa-head  | ghbpw-head-cn     |      2142.76  |          95.9596 |            0         |            2.92683   |
| GHB-Awanui-head   | ghbaw-head-hg     |      6226.52  |          95.82   |            0.194423  |            0.188602  |
| GHB-Spring-head   | ghbspring-head-hg |      1537.21  |          93.228  |            0         |            0.97561   |
| GHB-Poukawa-head  | ghbpw-head-hg     |       804.733 |          84.1078 |            0.0195122 |            0.35122   |
| GHB-Awanui-head   | ghbaw-head-cn     |       513.778 |          76.6077 |            0         |            0.487805  |
| GHB-Confined-head | ghbconf-head-hg   |       901.672 |          61.902  |            0.729136  |            0.0284524 |

### GHB Conductance

| Subcategory       | Pname             |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:------------------|:------------------|--------------:|-----------------:|---------------------:|---------------------:|
| GHB-Confined-cond | ghbconf-cond-cn   |    14354.8    |          99.3996 |             0        |             0.487805 |
| GHB-Confined-cond | ghbconf-cond-gr   |     3014.59   |          96.4425 |             0.574191 |             0.001714 |
| GHB-Spring-cond   | ghbspring-cond-cn |     2209.71   |          96.1574 |             0.487805 |             0        |
| GHB-Poukawa-cond  | ghbpw-cond-cn     |     1062.3    |          91.6977 |             1.46341  |             0        |
| GHB-Spring-cond   | ghbspring-cond-gr |     3568.42   |          88.2558 |             0.731707 |             0        |
| GHB-Awanui-cond   | ghbaw-cond-gr     |     2158.92   |          83.0095 |             0.115257 |             0.200244 |
| GHB-Poukawa-cond  | ghbpw-cond-gr     |     1523.08   |          14.0474 |             0.097561 |             0.097561 |
| GHB-Awanui-cond   | ghbaw-cond-cn     |       46.4548 |         -74.6444 |             0        |             0        |

### Hydraulic Conductivity

| Subcategory   | Pname         |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:--------------|:--------------|--------------:|-----------------:|---------------------:|---------------------:|
| K-pilot       | npfklayer1-pp |      9948.89  |          98.4773 |            0         |             5.29334  |
| K-grid        | npfklayer1-gr |      1857.17  |          74.6726 |            0.0051079 |             1.7137   |
| K-constant    | npfklayer1-cn |       204.305 |          24.1776 |            0         |             0.487805 |

### Storage

| Subcategory   | Pname          |   Mean_Change |   CV_Reduction_% |   Pct_at_Lower_Bound |   Pct_at_Upper_Bound |
|:--------------|:---------------|--------------:|-----------------:|---------------------:|---------------------:|
| S-pilot       | stosslayer1-pp |      24945.5  |          99.3882 |           0.00659196 |             0.494397 |
| S-constant    | stosslayer1-cn |       1673.5  |          94.7072 |           0          |             0        |
| S-grid        | stosslayer1-gr |       4465.26 |          92.2797 |           0.158345   |             0.117482 |

## Boundary Hitting Diagnostics

Parameters at bounds may indicate over-parameterization or improper bounds.

### Groups with >5% at Upper Bound:

| Subcategory   | Pname         |   Pct_at_Upper_Bound |
|:--------------|:--------------|---------------------:|
| K-pilot       | npfklayer1-pp |              5.29334 |

## Uncertainty Reduction Summary

| Category               |    mean |      min |     max |
|:-----------------------|--------:|---------:|--------:|
| GHB Conductance        | 61.7957 | -74.6444 | 99.3996 |
| GHB Head               | 88.3623 |  61.902  | 99.6891 |
| Hydraulic Conductivity | 65.7758 |  24.1776 | 98.4773 |
| Storage                | 95.4584 |  92.2797 | 99.3882 |

