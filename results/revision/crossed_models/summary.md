# Crossed-dependence reanalysis of real versus simulated performance

## F1

Units: 549 observations, 55 tools, 20 biological datasets, 11 benchmarks (6 with both real and simulated data); 5 values equal to 0 or 1.

- M0 submitted (benchmark fixed, tool random) [metric]: -0.334 (95% CI -0.370 to -0.299; Wald); p = 3.32e-76
- M1 crossed logit-normal (primary) [logit]: -1.996 (95% CI -2.567 to -1.166; benchmark cluster bootstrap (B=1999)); p = 1.54e-04
- M1 crossed logit-normal (primary) [logit]: -1.996 (95% CI -3.030 to -0.962; Wald (model-based)); p = 1.54e-04
- M1 crossed logit-normal (primary) [metric difference]: -0.460 (95% CI -0.539 to -0.278; benchmark cluster bootstrap (B=1999)); p = NA
- M2 crossed linear (sensitivity) [metric]: -0.352 (95% CI -0.422 to -0.218; benchmark cluster bootstrap (B=1999)); p = 1.31e-06
- M2 crossed linear (sensitivity) [metric]: -0.352 (95% CI -0.494 to -0.209; Wald (model-based)); p = 1.31e-06

- Within-benchmark contrast: mean real - simulated = -0.279 (median -0.331), negative in 5/6 benchmarks; Wilcoxon p = 0.062, sign test p = 0.219.
- Same tool within benchmark: mean paired difference -0.286, negative in 5/6 benchmarks.

- Leave-one-benchmark-out (M1, metric scale): -0.507 to -0.410; (M2): -0.388 to -0.310.
- Excluding Edgren-derived observations: M1 difference -0.489 (logit p = 1.39e-05); M2 -0.387 (95% CI -0.484 to -0.290; p = 5.23e-15).

## Precision

Units: 551 observations, 56 tools, 23 biological datasets, 11 benchmarks (7 with both real and simulated data); 61 values equal to 0 or 1.

- M0 submitted (benchmark fixed, tool random) [metric]: -0.393 (95% CI -0.429 to -0.357; Wald); p = 1.87e-100
- M1 crossed logit-normal (primary) [logit]: -3.161 (95% CI -4.104 to -2.090; benchmark cluster bootstrap (B=1999)); p = 1.43e-04
- M1 crossed logit-normal (primary) [logit]: -3.161 (95% CI -4.791 to -1.532; Wald (model-based)); p = 1.43e-04
- M1 crossed logit-normal (primary) [metric difference]: -0.474 (95% CI -0.724 to -0.245; benchmark cluster bootstrap (B=1999)); p = NA
- M2 crossed linear (sensitivity) [metric]: -0.403 (95% CI -0.525 to -0.258; benchmark cluster bootstrap (B=1999)); p = 2.99e-05
- M2 crossed linear (sensitivity) [metric]: -0.403 (95% CI -0.593 to -0.214; Wald (model-based)); p = 2.99e-05

- Within-benchmark contrast: mean real - simulated = -0.387 (median -0.397), negative in 6/7 benchmarks; Wilcoxon p = 0.031, sign test p = 0.125.
- Same tool within benchmark: mean paired difference -0.368, negative in 6/7 benchmarks.

- Leave-one-benchmark-out (M1, metric scale): -0.595 to -0.404; (M2): -0.479 to -0.341.
- Excluding Edgren-derived observations: M1 difference -0.492 (logit p = 3.99e-06); M2 -0.461 (95% CI -0.576 to -0.346; p = 4.32e-15).

## Recall

Units: 542 observations, 57 tools, 22 biological datasets, 11 benchmarks (8 with both real and simulated data); 22 values equal to 0 or 1.

- M0 submitted (benchmark fixed, tool random) [metric]: -0.123 (95% CI -0.164 to -0.082; Wald); p = 3.51e-09
- M1 crossed logit-normal (primary) [logit]: -0.338 (95% CI -0.829 to 0.098; benchmark cluster bootstrap (B=1999)); p = 0.485
- M1 crossed logit-normal (primary) [logit]: -0.338 (95% CI -1.287 to 0.610; Wald (model-based)); p = 0.485
- M1 crossed logit-normal (primary) [metric difference]: -0.079 (95% CI -0.198 to 0.020; benchmark cluster bootstrap (B=1999)); p = NA
- M2 crossed linear (sensitivity) [metric]: -0.095 (95% CI -0.154 to -0.021; benchmark cluster bootstrap (B=1999)); p = 0.111
- M2 crossed linear (sensitivity) [metric]: -0.095 (95% CI -0.213 to 0.022; Wald (model-based)); p = 0.111

- Within-benchmark contrast: mean real - simulated = -0.114 (median -0.084), negative in 7/8 benchmarks; Wilcoxon p = 0.055, sign test p = 0.070.
- Same tool within benchmark: mean paired difference -0.136, negative in 7/8 benchmarks.

- Leave-one-benchmark-out (M1, metric scale): -0.128 to -0.058; (M2): -0.118 to -0.060.
- Excluding Edgren-derived observations: M1 difference -0.060 (logit p = 0.644); M2 -0.095 (95% CI -0.228 to 0.038; p = 0.160).

All estimates are associations between benchmark settings and reported performance, not causal effects.
