# Precision-recall trade-off analysis

Tool-level means were calculated from observations with both recall/sensitivity and precision available.

## Correlations
- all: n=471, Pearson r=0.418 (p=2.31e-21), Spearman rho=0.333 (p=1.23e-13).
- without_edgren: n=391, Pearson r=0.449 (p=8.19e-21), Spearman rho=0.308 (p=4.72e-10).
- real: n=179, Pearson r=0.109 (p=1.46e-01), Spearman rho=0.094 (p=2.09e-01).
- real_without_edgren: n=99, Pearson r=0.118 (p=2.44e-01), Spearman rho=0.089 (p=3.81e-01).
- simulated: n=264, Pearson r=0.645 (p=1.93e-32), Spearman rho=0.303 (p=5.44e-07).

## High recall / low precision profile
- No tools met the quartile-based high recall / low precision definition.

## Low recall / high precision profile
- JAFFA-Assembly: mean recall=0.555, mean precision=0.807, mean F1=0.642, n=5.
- TrinityFusion-D: mean recall=0.412, mean precision=0.961, mean F1=0.575, n=3.

## Largest recall-minus-precision gaps
- STAR-SEQR: recall - precision = 0.258; mean recall=0.842; mean precision=0.583; n=11.
- FuSeq: recall - precision = 0.231; mean recall=0.685; mean precision=0.454; n=10.
- JAFFA: recall - precision = 0.212; mean recall=0.567; mean precision=0.355; n=13.
- GFusion: recall - precision = 0.210; mean recall=0.779; mean precision=0.569; n=3.
- EricScript: recall - precision = 0.192; mean recall=0.602; mean precision=0.410; n=18.
- ChimeraScan: recall - precision = 0.183; mean recall=0.578; mean precision=0.394; n=29.
- MetaFusion: recall - precision = 0.162; mean recall=0.920; mean precision=0.757; n=6.
- Arriba: recall - precision = 0.146; mean recall=0.797; mean precision=0.652; n=23.
- ChimeRScope: recall - precision = 0.131; mean recall=0.254; mean precision=0.123; n=7.
- Pizzly: recall - precision = 0.111; mean recall=0.663; mean precision=0.552; n=10.
