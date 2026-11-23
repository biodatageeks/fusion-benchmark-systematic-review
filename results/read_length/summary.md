# Read-length analysis

Input observations with 50/75/100 bp and F1: 392

## Mean F1 by read length
- 50 bp: n=160, mean=0.514, median=0.537.
- 75 bp: n=107, mean=0.677, median=0.788.
- 100 bp: n=125, mean=0.607, median=0.716.

Kruskal-Wallis H=27.86, p=8.93e-07.

## Mixed-effects model with crossed random intercepts (tool, biological dataset group, benchmark), reference 100 bp
n=392 observations, 38 tools, 16 dataset groups, 12 benchmarks; converged: M1 True, M2 True.
- 75 bp: M1 difference in F1=0.055, 95% CI -0.094 to 0.189, p=0.4637; M2 difference=0.030, 95% CI -0.052 to 0.112, p=0.4723.
- 50 bp: M1 difference in F1=-0.063, 95% CI -0.181 to 0.058, p=0.3125; M2 difference=-0.029, 95% CI -0.094 to 0.036, p=0.3801.
