# New-tool performance within original benchmarks

Benchmark-level comparison based on mean tool performance within each benchmark.

New tool/family ranked first by mean F1 in 6 of 7 new-tool benchmarks with available F1.
Using F1 where available and recall otherwise, the new tool/family ranked first in 8 of 9 new-tool benchmarks.
New tool/family had higher mean F1 than the average of comparator tools in 7 of 7 benchmarks with available F1.

| Benchmark | Study | New tool | best new tool | primary metric/rank | best overall | new mean F1 | others mean F1 | datasets best |
|---:|---|---|---|---:|---|---:|---:|---:|
| 1 | Haas 2019 | STAR-Fusion; TrinityFusion | STAR-Fusion | F1 1/24 | STAR-Fusion | 0.743 | 0.717 | 1/3 |
| 2 | Uhrig 2021 | Arriba | Arriba | F1 1/7 | Arriba | 0.802 | 0.469 | 14/29 |
| 4 | Balan 2021 | SeekFusion | SeekFusion | Recall_Sensitivity 1/5 | SeekFusion | nan | nan | NA |
| 5 | Zhao 2017 | GFusion | GFusion | F1 1/5 | GFusion | 0.682 | 0.494 | 3/5 |
| 6 | Vu 2018 | FuSeq | FuSeq | F1 2/7 | FuSeq+assembly | 0.454 | 0.353 | 4/9 |
| 8 | Apostolides 2021 | MetaFusion | MetaFusion-cf4 | F1 1/14 | MetaFusion-cf4 | 0.808 | 0.593 | 6/6 |
| 9 | Zhang 2016 | INTEGRATE | INTEGRATE (combined; 1 read) | Recall_Sensitivity 1/13 | INTEGRATE (combined; 1 read) | nan | nan | NA |
| 10 | Davidson 2015 | JAFFA | JAFFA-Hybrid | F1 1/12 | JAFFA-Hybrid | 0.384 | 0.359 | 3/10 |
| 11 | Vicente-Garces 2023 | Fusion InPipe | Fusion InPipe-3 | F1 1/8 | Fusion InPipe-3 | 0.833 | 0.823 | 2/2 |

## Notes
- Benchmark 4: recall only (SeekFusion)
- Benchmark 9: recall only (INTEGRATE)
