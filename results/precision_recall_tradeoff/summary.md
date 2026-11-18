# Precision-recall trade-off analysis

Tool-level means were calculated from observations with both recall/sensitivity and precision available.

## Correlations
- all: n=662, Pearson r=0.164 (p=2.23e-05), Spearman rho=0.089 (p=2.25e-02).
- without_edgren: n=588, Pearson r=0.093 (p=2.48e-02), Spearman rho=0.012 (p=7.78e-01).
- real: n=227, Pearson r=0.010 (p=8.86e-01), Spearman rho=-0.024 (p=7.20e-01).
- real_without_edgren: n=180, Pearson r=-0.120 (p=1.09e-01), Spearman rho=-0.145 (p=5.14e-02).
- simulated: n=252, Pearson r=0.260 (p=2.99e-05), Spearman rho=0.038 (p=5.49e-01).

## High recall / low precision profile
- FuSeq: mean recall=0.748, mean precision=0.413, mean F1=0.495, n=15.

## Low recall / high precision profile
- PRADA: mean recall=0.475, mean precision=0.935, mean F1=0.592, n=30.
- FusionHunter: mean recall=0.395, mean precision=0.900, mean F1=0.543, n=4.
- TrinityFusion-D: mean recall=0.413, mean precision=0.960, mean F1=0.576, n=3.
- MetaFusion-cf7: mean recall=0.369, mean precision=1.000, mean F1=0.527, n=6.

## Largest recall-minus-precision gaps
- JAFFA-Direct (HC+MC): recall - precision = 0.496; mean recall=0.687; mean precision=0.191; n=7.
- Hera: recall - precision = 0.444; mean recall=0.457; mean precision=0.012; n=3.
- deFuse: recall - precision = 0.355; mean recall=0.725; mean precision=0.370; n=50.
- JAFFA-Direct: recall - precision = 0.349; mean recall=0.713; mean precision=0.363; n=10.
- FuSeq: recall - precision = 0.335; mean recall=0.748; mean precision=0.413; n=15.
- InFusion: recall - precision = 0.253; mean recall=0.600; mean precision=0.347; n=39.
- STAR-SEQR: recall - precision = 0.245; mean recall=0.843; mean precision=0.598; n=11.
- GFusion: recall - precision = 0.238; mean recall=0.799; mean precision=0.561; n=6.
- Arriba: recall - precision = 0.202; mean recall=0.889; mean precision=0.687; n=50.
- SOAPfuse: recall - precision = 0.190; mean recall=0.659; mean precision=0.468; n=57.
