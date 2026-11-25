# New-tool performance within original benchmarks

Benchmark-level comparison based on mean tool performance within each benchmark.

New tool/family ranked first by mean F1 in 12 of 12 new-tool benchmarks with available F1.
Using F1 where available and recall otherwise, the new tool/family ranked first in 13 of 13 new-tool benchmarks.
New tool/family had higher mean F1 than the average of comparator tools in 12 of 12 benchmarks with available F1.

| Benchmark | Study | New tool | best new tool | primary metric/rank | best overall | new mean F1 | others mean F1 | datasets best |
|---:|---|---|---|---:|---|---:|---:|---:|
| 1 | Haas 2019 | STAR-Fusion; TrinityFusion | STAR-Fusion | F1 1/25 | STAR-Fusion | 0.742 | 0.715 | 1/3 |
| 2 | Uhrig 2021 | Arriba | Arriba | F1 1/7 | Arriba | 0.864 | 0.363 | 2/2 |
| 4 | Balan 2021 | SeekFusion | SeekFusion | Recall_Sensitivity 1/4 | SeekFusion | nan | nan | NA |
| 5 | Zhao 2017 | GFusion | GFusion | F1 1/5 | GFusion | 0.825 | 0.475 | 1/2 |
| 6 | Vu 2018 | FuSeq | FuSeq | F1 1/6 | FuSeq | 0.432 | 0.257 | 2/4 |
| 8 | Apostolides 2021 | MetaFusion | MetaFusion.top_3 | F1 1/9 | MetaFusion.top_3 | 0.832 | 0.593 | 6/6 |
| 9 | Zhang 2016 | INTEGRATE | INTEGRATE C | F1 1/12 | INTEGRATE C | 0.694 | 0.206 | 1/1 |
| 10 | Davidson 2015 | JAFFA | JAFFA-Hybrid | F1 1/7 | JAFFA-Hybrid | 0.795 | 0.529 | 1/2 |
| 11 | Vicente-Garces 2023 | Fusion InPipe | Fusion InPipe-3 | F1 1/8 | Fusion InPipe-3 | 0.833 | 0.823 | 2/2 |
| 15 | Rodriguez-Martin 2017 | ChimPipe | ChimPipe | F1 1/5 | ChimPipe | 0.778 | 0.518 | 6/16 |
| 16 | Jia 2013 | SOAPfuse | SOAPfuse | F1 1/3 | SOAPfuse | 0.930 | 0.792 | 8/9 |
| 17 | Benelli 2012 | EricScript | EricScript | F1 1/5 | EricScript | 0.797 | 0.538 | 1/1 |
| 18 | Nicorici 2014 | FusionCatcher | FusionCatcher | F1 1/7 | FusionCatcher | 0.892 | 0.443 | 1/1 |

## Notes
- Benchmark 4: recall only (SeekFusion)
- Benchmark 14: new tool not found in observations (SimFuse)
