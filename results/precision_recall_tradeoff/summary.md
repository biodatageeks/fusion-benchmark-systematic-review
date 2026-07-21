# Precision-recall trade-off analysis

Tool-level means were calculated from observations with both recall/sensitivity and precision available.

## Correlations
- all: n=315, Pearson r=0.417 (p=1.15e-14), Spearman rho=0.355 (p=8.63e-11).
- without_edgren: n=250, Pearson r=0.441 (p=2.51e-13), Spearman rho=0.321 (p=2.09e-07).
- real: n=150, Pearson r=0.129 (p=1.16e-01), Spearman rho=0.169 (p=3.89e-02).
- real_without_edgren: n=85, Pearson r=0.123 (p=2.64e-01), Spearman rho=0.153 (p=1.64e-01).
- simulated: n=158, Pearson r=0.568 (p=7.11e-15), Spearman rho=0.330 (p=2.25e-05).

## High recall / low precision profile
- No tools met the quartile-based high recall / low precision definition.

## Low recall / high precision profile
- FusionHunter: mean recall=0.337, mean precision=0.860, mean F1=0.483, n=3.
- TrinityFusion-D: mean recall=0.427, mean precision=0.960, mean F1=0.590, n=3.

## Largest recall-minus-precision gaps
- ChimeraScan: recall - precision = 0.234; mean recall=0.672; mean precision=0.438; n=11.
- FuSeq: recall - precision = 0.221; mean recall=0.610; mean precision=0.389; n=11.
- JAFFA: recall - precision = 0.196; mean recall=0.526; mean precision=0.330; n=14.
- STAR-SEQR: recall - precision = 0.196; mean recall=0.880; mean precision=0.684; n=9.
- deFuse: recall - precision = 0.193; mean recall=0.679; mean precision=0.486; n=13.
- EricScript: recall - precision = 0.174; mean recall=0.579; mean precision=0.404; n=16.
- Arriba: recall - precision = 0.163; mean recall=0.790; mean precision=0.627; n=19.
- MetaFusion: recall - precision = 0.162; mean recall=0.920; mean precision=0.757; n=6.
- ChimeRScope: recall - precision = 0.130; mean recall=0.253; mean precision=0.123; n=7.
- FusionMap: recall - precision = 0.109; mean recall=0.593; mean precision=0.484; n=12.
