# Crossed-dependence models, beta regression, real-data subtypes

## F1

518 observations (216 real, 302 simulated), 49 tools, 22 biological dataset groups, 16 benchmarks (7 with both real and simulated data).

- M0 benchmark fixed, tool random (ignores dataset dependence) [metric difference]: -0.338 (95% CI -0.378 to -0.297; Wald); p = 1.9e-60
- M1 crossed logit-normal (primary) [logit]: -1.728 (95% CI -2.731 to -0.768; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M1 crossed logit-normal (primary) [logit]: -1.728 (95% CI -2.824 to -0.632; Wald); p = 0.002
- M1 crossed logit-normal (primary) [metric difference]: -0.407 (95% CI -0.566 to -0.179; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M1 crossed logit-normal (primary) [metric difference]: -0.407 (95% CI -0.581 to -0.145; Wald); p = 0.002
- M2 crossed linear (sensitivity) [metric difference]: -0.310 (95% CI -0.434 to -0.143; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M2 crossed linear (sensitivity) [metric difference]: -0.310 (95% CI -0.459 to -0.161; Wald); p = 4.7e-05
- M3 beta regression, tool fixed [logit]: -1.216 (95% CI -1.741 to -0.834; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M3 beta regression, tool fixed [logit]: -1.216 (95% CI -1.407 to -1.026; Wald); p = 4.8e-36
- M3 beta regression, tool fixed [metric difference]: -0.274 (95% CI -0.370 to -0.183; benchmark cluster bootstrap (B=999)); p = 0.0e+00

- Convergence (after optimizer fallback): bootstrap M1 999/999, M2 999/999 replicates; leave-one-benchmark-out M1 16/16 fits.
- Leave-one-benchmark-out, M1 metric difference: -0.514 to -0.364; M2: -0.374 to -0.253.
- Without Edgren-derived observations: M1 -0.405 (logit p = 0.0082); M2 -0.391 (95% CI -0.480 to -0.303; p = 4.9e-18).
- Within-benchmark contrast (real − simulated, one value per benchmark): mean -0.281, negative in 6/7 benchmarks (sign test p = 0.12); same tool within benchmark: mean -0.285, negative in 6/7.

- cell line vs simulated (n = 153, 15 datasets, 11 benchmarks): M1 -0.310 (95% CI -0.560 to 0.013; p = 0.062); M2 -0.255 (p = 0.0066)
- patient-derived vs simulated (n = 25, 3 datasets, 3 benchmarks): M1 -0.503 (95% CI -0.658 to -0.124; p = 0.0099); M2 -0.332 (p = 0.0056)
- mixed vs simulated (n = 38, 5 datasets, 4 benchmarks): M1 -0.499 (95% CI -0.681 to 0.063; p = 0.087); M2 -0.389 (p = 0.017)
- patient-derived vs cell line: M1 logit -0.958 (p = 0.32); M2 -0.077 (95% CI -0.320 to 0.166; p = 0.54)

## Precision

452 observations (188 real, 264 simulated), 47 tools, 22 biological dataset groups, 16 benchmarks (8 with both real and simulated data).

- M0 benchmark fixed, tool random (ignores dataset dependence) [metric difference]: -0.420 (95% CI -0.466 to -0.374; Wald); p = 9.1e-72
- M1 crossed logit-normal (primary) [logit]: -2.906 (95% CI -4.316 to -1.493; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M1 crossed logit-normal (primary) [logit]: -2.906 (95% CI -4.531 to -1.281; Wald); p = 4.6e-04
- M1 crossed logit-normal (primary) [metric difference]: -0.520 (95% CI -0.774 to -0.243; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M1 crossed logit-normal (primary) [metric difference]: -0.520 (95% CI -0.807 to -0.149; Wald); p = 4.6e-04
- M2 crossed linear (sensitivity) [metric difference]: -0.384 (95% CI -0.565 to -0.246; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M2 crossed linear (sensitivity) [metric difference]: -0.384 (95% CI -0.593 to -0.175; Wald); p = 3.1e-04
- M3 beta regression, tool fixed [logit]: -1.543 (95% CI -2.103 to -1.168; benchmark cluster bootstrap (B=999)); p = 0.0e+00
- M3 beta regression, tool fixed [logit]: -1.543 (95% CI -1.777 to -1.310; Wald); p = 2.6e-38
- M3 beta regression, tool fixed [metric difference]: -0.327 (95% CI -0.420 to -0.240; benchmark cluster bootstrap (B=999)); p = 0.0e+00

- Convergence (after optimizer fallback): bootstrap M1 999/999, M2 999/999 replicates; leave-one-benchmark-out M1 16/16 fits.
- Leave-one-benchmark-out, M1 metric difference: -0.639 to -0.459; M2: -0.460 to -0.355.
- Without Edgren-derived observations: M1 -0.553 (logit p = 1.0e-07); M2 -0.438 (95% CI -0.581 to -0.296; p = 1.5e-09).
- Within-benchmark contrast (real − simulated, one value per benchmark): mean -0.417, negative in 7/8 benchmarks (sign test p = 0.07); same tool within benchmark: mean -0.394, negative in 7/8.

- cell line vs simulated (n = 149, 17 datasets, 12 benchmarks): M1 -0.486 (95% CI -0.824 to -0.074; p = 0.0069); M2 -0.350 (p = 0.0053)
- patient-derived vs simulated (n = 15, 3 datasets, 3 benchmarks): M1 -0.472 (95% CI -0.877 to 0.013; p = 0.068); M2 -0.398 (p = 0.028)
- mixed vs simulated (n = 24, 4 datasets, 3 benchmarks): M1 -0.703 (95% CI -0.916 to 0.011; p = 0.061); M2 -0.538 (p = 0.038)
- patient-derived vs cell line: M1 logit 0.055 (p = 0.97); M2 -0.048 (95% CI -0.414 to 0.318; p = 0.8)

## Recall

485 observations (209 real, 276 simulated), 49 tools, 23 biological dataset groups, 17 benchmarks (9 with both real and simulated data).

- M0 benchmark fixed, tool random (ignores dataset dependence) [metric difference]: -0.140 (95% CI -0.189 to -0.092; Wald); p = 1.8e-08
- M1 crossed logit-normal (primary) [logit]: -0.443 (95% CI -1.238 to 0.134; benchmark cluster bootstrap (B=999)); p = 0.13
- M1 crossed logit-normal (primary) [logit]: -0.443 (95% CI -1.509 to 0.624; Wald); p = 0.42
- M1 crossed logit-normal (primary) [metric difference]: -0.102 (95% CI -0.298 to 0.029; benchmark cluster bootstrap (B=999)); p = 0.13
- M1 crossed logit-normal (primary) [metric difference]: -0.102 (95% CI -0.360 to 0.118; Wald); p = 0.42
- M2 crossed linear (sensitivity) [metric difference]: -0.058 (95% CI -0.190 to 0.018; benchmark cluster bootstrap (B=999)); p = 0.13
- M2 crossed linear (sensitivity) [metric difference]: -0.058 (95% CI -0.210 to 0.093; Wald); p = 0.45
- M3 beta regression, tool fixed [logit]: -0.416 (95% CI -0.907 to 0.089; benchmark cluster bootstrap (B=999)); p = 0.094
- M3 beta regression, tool fixed [logit]: -0.416 (95% CI -0.617 to -0.214; Wald); p = 5.3e-05
- M3 beta regression, tool fixed [metric difference]: -0.096 (95% CI -0.201 to 0.020; benchmark cluster bootstrap (B=999)); p = 0.094

- Convergence (after optimizer fallback): bootstrap M1 999/999, M2 999/999 replicates; leave-one-benchmark-out M1 17/17 fits.
- Leave-one-benchmark-out, M1 metric difference: -0.157 to -0.073; M2: -0.113 to -0.036.
- Without Edgren-derived observations: M1 -0.093 (logit p = 0.51); M2 -0.056 (95% CI -0.222 to 0.110; p = 0.51).
- Within-benchmark contrast (real − simulated, one value per benchmark): mean -0.110, negative in 6/9 benchmarks (sign test p = 0.51); same tool within benchmark: mean -0.108, negative in 6/9.

- cell line vs simulated (n = 159, 17 datasets, 13 benchmarks): M1 -0.180 (95% CI -0.459 to 0.099; p = 0.24); M2 -0.102 (p = 0.27)
- patient-derived vs simulated (n = 22, 4 datasets, 4 benchmarks): M1 0.051 (95% CI -0.340 to 0.255; p = 0.77); M2 0.047 (p = 0.69)
- mixed vs simulated (n = 28, 5 datasets, 4 benchmarks): M1 -0.080 (95% CI -0.510 to 0.236; p = 0.73); M2 -0.071 (p = 0.63)
- patient-derived vs cell line: M1 logit 1.002 (p = 0.25); M2 0.148 (95% CI -0.094 to 0.391; p = 0.23)

## Benchmark purpose (descriptive)

- F1, independent_comparison: 209 observations, 4 benchmarks, median 0.472 (IQR 0.209–0.716), mean 0.455
- F1, new_tool: 309 observations, 12 benchmarks, median 0.722 (IQR 0.411–0.850), mean 0.623
- Precision, independent_comparison: 134 observations, 4 benchmarks, median 0.530 (IQR 0.260–0.810), mean 0.516
- Precision, new_tool: 318 observations, 12 benchmarks, median 0.861 (IQR 0.438–0.970), mean 0.690
- Recall, independent_comparison: 158 observations, 4 benchmarks, median 0.597 (IQR 0.264–0.764), mean 0.507
- Recall, new_tool: 327 observations, 13 benchmarks, median 0.710 (IQR 0.524–0.832), mean 0.648
- Same tool, F1 (all tools): higher in new-tool benchmarks for 18/21 tools; median difference 0.172
- Same tool, F1 (without the authors' own new tool): higher in new-tool benchmarks for 17/20 tools; median difference 0.133
- Same tool, Precision (all tools): higher in new-tool benchmarks for 13/16 tools; median difference 0.142
- Same tool, Precision (without the authors' own new tool): higher in new-tool benchmarks for 12/15 tools; median difference 0.130
- Same tool, Recall (all tools): higher in new-tool benchmarks for 15/19 tools; median difference 0.147
- Same tool, Recall (without the authors' own new tool): higher in new-tool benchmarks for 13/18 tools; median difference 0.126

All estimates are associations between benchmark settings and reported performance, not causal effects.
