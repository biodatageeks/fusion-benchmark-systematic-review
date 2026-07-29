# Precision-recall trade-off analysis

Tool-level means were calculated from observations with both recall/sensitivity and precision available.

## Correlations
- all: n=324, Pearson r=0.335 (p=6.18e-10), Spearman rho=0.290 (p=1.07e-07).
- without_edgren: n=257, Pearson r=0.352 (p=6.27e-09), Spearman rho=0.248 (p=5.81e-05).
- real: n=164, Pearson r=0.126 (p=1.08e-01), Spearman rho=0.127 (p=1.06e-01).
- real_without_edgren: n=97, Pearson r=0.163 (p=1.10e-01), Spearman rho=0.135 (p=1.86e-01).
- simulated: n=154, Pearson r=0.445 (p=7.17e-09), Spearman rho=0.278 (p=4.80e-04).

## High recall / low precision profile
- No tools met the quartile-based high recall / low precision definition.

## Low recall / high precision profile
- PRADA: mean recall=0.535, mean precision=0.967, mean F1=0.684, n=3.
- FusionHunter: mean recall=0.336, mean precision=0.859, mean F1=0.481, n=3.
- TrinityFusion-D: mean recall=0.427, mean precision=0.960, mean F1=0.590, n=3.

## Largest recall-minus-precision gaps
- FuSeq: recall - precision = 0.244; mean recall=0.671; mean precision=0.428; n=10.
- ChimeraScan: recall - precision = 0.233; mean recall=0.672; mean precision=0.438; n=11.
- JAFFA: recall - precision = 0.212; mean recall=0.567; mean precision=0.355; n=13.
- STAR-SEQR: recall - precision = 0.196; mean recall=0.880; mean precision=0.684; n=9.
- EricScript: recall - precision = 0.175; mean recall=0.579; mean precision=0.404; n=16.
- MetaFusion: recall - precision = 0.162; mean recall=0.920; mean precision=0.757; n=6.
- deFuse: recall - precision = 0.139; mean recall=0.698; mean precision=0.559; n=15.
- Arriba: recall - precision = 0.131; mean recall=0.791; mean precision=0.660; n=21.
- ChimeRScope: recall - precision = 0.131; mean recall=0.254; mean precision=0.123; n=7.
- FusionMap: recall - precision = 0.111; mean recall=0.594; mean precision=0.483; n=12.
