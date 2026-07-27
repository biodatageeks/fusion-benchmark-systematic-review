# New-tool performance within original benchmarks

Benchmark-level comparison based on mean tool performance within each benchmark.

New tool/family ranked first by mean F1 in 7 of 7 new-tool benchmarks with available F1.
Using F1 where available and recall otherwise, the new tool/family ranked first in 8 of 8 new-tool benchmarks.
New tool/family had higher mean F1 than the average of comparator tools in 7 of 7 benchmarks with available F1.

| Benchmark | Study | New tool | best new tool | primary metric/rank | best overall | new mean F1 | others mean F1 | datasets best |
|---:|---|---|---|---:|---|---:|---:|---:|
| 1 | Haas 2019 | STAR-Fusion; TrinityFusion | STAR-Fusion | F1 1/25 | STAR-Fusion | 0.765 | 0.727 | 1/3 |
| 2 | Uhrig 2021 | Arriba | Arriba | F1 1/7 | Arriba | 0.684 | 0.318 | 3/3 |
| 4 | Balan 2021 | SeekFusion | SeekFusion | Recall_Sensitivity 1/4 | SeekFusion | nan | nan | NA |
| 5 | Zhao 2017 | GFusion | GFusion | F1 1/5 | GFusion | 0.825 | 0.475 | 1/2 |
| 6 | Vu 2018 | FuSeq | FuSeq | F1 1/6 | FuSeq | 0.385 | 0.258 | 2/4 |
| 8 | Apostolides 2021 | MetaFusion | MetaFusion.top_3 | F1 1/9 | MetaFusion.top_3 | 0.832 | 0.593 | 6/6 |
| 9 | Zhang 2016 | INTEGRATE | INTEGRATE C | F1 1/12 | INTEGRATE C | 0.694 | 0.262 | 1/1 |
| 10 | Davidson 2015 | JAFFA | JAFFA-Hybrid | F1 1/7 | JAFFA-Hybrid | 0.795 | 0.529 | 1/2 |

## Notes
- Benchmark 4: recall only (SeekFusion)
