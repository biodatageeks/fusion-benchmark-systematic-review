"""Real-versus-simulated analysis respecting crossed dependence (R1.1, Editor E.1).

Each tool score is an observation of (tool x biological dataset x benchmark publication).
Tools recur across benchmarks and some biological datasets (above all Edgren) recur across
publications, so these three factors are crossed rather than nested.

Models (fitted separately for F1, precision and recall):

  M1  primary   logit(metric*) ~ real + (1|tool) + (1|bio_dataset) + (1|benchmark)
                metric* = Smithson-Verkuilen squeezed metric, so 0/1 values are allowed.
                Effect reported on the logit scale and as a difference on the metric scale
                (expit(b0 + b1) - expit(b0), i.e. for a typical tool/dataset/benchmark).
  M2  sens.     metric ~ real + (1|tool) + (1|bio_dataset) + (1|benchmark)   (identity scale)
  M0  original  metric ~ real + C(benchmark) + (1|tool)  (submitted model, for comparison only)
  W   benchmark-level within-benchmark contrasts: mean(real) - mean(simulated) per benchmark,
                one value per benchmark (no dependence on the number of tool scores).

Uncertainty for M1/M2: model-based Wald intervals and a cluster bootstrap that resamples
whole benchmarks with replacement (primary CI). The number of clusters is reported with
every estimate.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from revision_utils import (
    METRICS,
    REVISION_DIR,
    expit,
    fmt,
    fmt_p,
    load_observations,
    logit,
    squeeze,
    unit_counts,
    write_latex_table,
)

OUT = REVISION_DIR / "crossed_models"
VC = {"tool": "0 + C(Tool)", "bio": "0 + C(bio_dataset_id)", "bench": "0 + C(Benchmark)"}


def prepare(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    d = df.dropna(subset=[metric]).copy()
    d["y"] = d[metric].astype(float)
    d["y_logit"] = logit(squeeze(d["y"]))
    d["one"] = 1
    return d


def fit_crossed(d: pd.DataFrame, response: str):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return smf.mixedlm(f"{response} ~ real", d, groups=d["one"], vc_formula=VC).fit(reml=True)


def crossed_estimates(d: pd.DataFrame) -> dict[str, float]:
    m1 = fit_crossed(d, "y_logit")
    m2 = fit_crossed(d, "y")
    b0, b1 = m1.params["Intercept"], m1.params["real"]
    return {
        "m1_logit": float(b1),
        "m1_diff": float(expit(b0 + b1) - expit(b0)),
        "m2_diff": float(m2.params["real"]),
    }


def _bootstrap_worker(args: tuple[pd.DataFrame, int]) -> dict[str, float] | None:
    d, seed = args
    rng = np.random.default_rng(seed)
    benches = d["Benchmark"].unique()
    pick = rng.choice(benches, size=len(benches), replace=True)
    parts = [d[d["Benchmark"] == b].assign(Benchmark=f"{b}#{i}") for i, b in enumerate(pick)]
    s = pd.concat(parts, ignore_index=True)
    if s["real"].nunique() < 2:
        return None
    try:
        return crossed_estimates(s)
    except Exception:
        return None


def cluster_bootstrap(d: pd.DataFrame, n_boot: int, seed: int, workers: int) -> pd.DataFrame:
    jobs = [(d, seed + i) for i in range(n_boot)]
    if workers > 1:
        with mp.get_context("fork").Pool(workers) as pool:
            results = pool.map(_bootstrap_worker, jobs, chunksize=8)
    else:
        results = [_bootstrap_worker(job) for job in jobs]
    return pd.DataFrame([r for r in results if r is not None])


def within_benchmark(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    d = df.dropna(subset=[metric])
    rows = []
    for bench, g in d.groupby("Benchmark"):
        if g["real"].nunique() < 2:
            continue
        sim, real = g.loc[g["real"] == 0, metric], g.loc[g["real"] == 1, metric]
        tools_both = sorted(set(g.loc[g["real"] == 0, "Tool"]) & set(g.loc[g["real"] == 1, "Tool"]))
        paired = [
            g.loc[(g["Tool"] == t) & (g["real"] == 1), metric].mean() - g.loc[(g["Tool"] == t) & (g["real"] == 0), metric].mean()
            for t in tools_both
        ]
        rows.append({
            "metric": metric,
            "benchmark": bench,
            "benchmark_label": g["Benchmark_label"].iloc[0],
            "n_sim": len(sim),
            "n_real": len(real),
            "mean_sim": sim.mean(),
            "mean_real": real.mean(),
            "diff_real_minus_sim": real.mean() - sim.mean(),
            "n_tools_both": len(tools_both),
            "paired_tool_diff_mean": float(np.mean(paired)) if paired else np.nan,
        })
    return pd.DataFrame(rows)


def summarize_within(w: pd.DataFrame) -> dict[str, object]:
    diffs = w["diff_real_minus_sim"].to_numpy()
    k = len(diffs)
    out = {"k_benchmarks": k, "mean_diff": float(np.mean(diffs)) if k else np.nan,
           "median_diff": float(np.median(diffs)) if k else np.nan, "n_negative": int((diffs < 0).sum())}
    out["wilcoxon_p"] = float(stats.wilcoxon(diffs).pvalue) if k >= 5 else np.nan
    out["sign_test_p"] = float(stats.binomtest(int((diffs < 0).sum()), k, 0.5).pvalue) if k else np.nan
    paired = w["paired_tool_diff_mean"].dropna().to_numpy()
    out["k_paired"] = len(paired)
    out["mean_paired_tool_diff"] = float(np.mean(paired)) if len(paired) else np.nan
    out["n_negative_paired"] = int((paired < 0).sum())
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--n-boot", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20261117)
    parser.add_argument("--workers", type=int, default=max(1, mp.cpu_count() - 1))
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    df = load_observations()
    model_rows, lobo_rows, within_rows, within_summary, vc_rows = [], [], [], [], []

    for metric in METRICS:
        d = prepare(df, metric)
        counts = unit_counts(d)
        k_both = int(d.groupby("Benchmark")["real"].nunique().eq(2).sum())
        n_boundary = int(((d["y"] == 0) | (d["y"] == 1)).sum())

        # M0: submitted specification, reproduced for comparison
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m0 = smf.mixedlm("y ~ real + C(Benchmark)", d, groups=d["Tool"]).fit(reml=True)
        ci0 = m0.conf_int().loc["real"]
        model_rows.append({"metric": metric, "model": "M0 submitted (benchmark fixed, tool random)", "scale": "metric",
                           "estimate": m0.params["real"], "ci_low": ci0[0], "ci_high": ci0[1], "ci_method": "Wald",
                           "p": m0.pvalues["real"], **counts, "k_benchmarks_both_types": k_both, "n_boundary_values": n_boundary})

        m1 = fit_crossed(d, "y_logit")
        m2 = fit_crossed(d, "y")
        est = crossed_estimates(d)
        boot = cluster_bootstrap(d, args.n_boot, args.seed, args.workers)
        boot.to_csv(OUT / f"bootstrap_{metric.lower()}.csv", index=False)
        q = lambda col: np.nanpercentile(boot[col], [2.5, 97.5]) if len(boot) else (np.nan, np.nan)

        ci1 = m1.conf_int().loc["real"]
        for name, scale, value, wald, p, col in [
            ("M1 crossed logit-normal (primary)", "logit", est["m1_logit"], ci1, m1.pvalues["real"], "m1_logit"),
            ("M1 crossed logit-normal (primary)", "metric difference", est["m1_diff"], (np.nan, np.nan), np.nan, "m1_diff"),
            ("M2 crossed linear (sensitivity)", "metric", est["m2_diff"], m2.conf_int().loc["real"], m2.pvalues["real"], "m2_diff"),
        ]:
            lo, hi = q(col)
            base = {"metric": metric, "model": name, "scale": scale, "estimate": value, "p": p, **counts,
                    "k_benchmarks_both_types": k_both, "n_boundary_values": n_boundary}
            model_rows.append({**base, "ci_low": lo, "ci_high": hi, "ci_method": f"benchmark cluster bootstrap (B={len(boot)})",
                               "p_boot_two_sided": float(2 * min((boot[col] >= 0).mean(), (boot[col] <= 0).mean())) if len(boot) else np.nan})
            if np.isfinite(wald[0]):
                model_rows.append({**base, "ci_low": wald[0], "ci_high": wald[1], "ci_method": "Wald (model-based)"})

        vc_rows.append({"metric": metric, "model": "M1", **dict(zip(m1.model.exog_vc.names, m1.vcomp)), "residual": m1.scale})
        vc_rows.append({"metric": metric, "model": "M2", **dict(zip(m2.model.exog_vc.names, m2.vcomp)), "residual": m2.scale})

        # leave-one-benchmark-out and Edgren exclusion
        for bench in sorted(d["Benchmark"].unique(), key=lambda x: int(x)):
            s = d[d["Benchmark"] != bench]
            if s["real"].nunique() < 2:
                continue
            e = crossed_estimates(s)
            lobo_rows.append({"metric": metric, "excluded": bench, **e})
        s = d[d["Is_Edgren_Binary"] != "yes"]
        e = crossed_estimates(s)
        m1e, m2e = fit_crossed(s, "y_logit"), fit_crossed(s, "y")
        lobo_rows.append({"metric": metric, "excluded": "Edgren-derived observations", **e,
                          "m1_logit_p": m1e.pvalues["real"], "m2_p": m2e.pvalues["real"],
                          "m2_ci_low": m2e.conf_int().loc["real"][0], "m2_ci_high": m2e.conf_int().loc["real"][1],
                          **{f"noedgren_{k}": v for k, v in unit_counts(s).items()}})

        w = within_benchmark(df, metric)
        within_rows.append(w)
        within_summary.append({"metric": metric, **summarize_within(w)})

    models = pd.DataFrame(model_rows)
    lobo = pd.DataFrame(lobo_rows)
    within = pd.concat(within_rows, ignore_index=True)
    wsum = pd.DataFrame(within_summary)
    models.to_csv(OUT / "models.csv", index=False)
    lobo.to_csv(OUT / "sensitivity_lobo_edgren.csv", index=False)
    within.to_csv(OUT / "within_benchmark_contrasts.csv", index=False)
    wsum.to_csv(OUT / "within_benchmark_summary.csv", index=False)
    pd.DataFrame(vc_rows).to_csv(OUT / "variance_components.csv", index=False)

    # Markdown summary
    lines = ["# Crossed-dependence reanalysis of real versus simulated performance", ""]
    for metric in METRICS:
        mm = models[models.metric == metric]
        c = mm.iloc[0]
        lines += [f"## {metric}", "",
                  f"Units: {c.n_obs} observations, {c.n_tools} tools, {c.n_bio_datasets} biological datasets, "
                  f"{c.n_benchmarks} benchmarks ({c.k_benchmarks_both_types} with both real and simulated data); "
                  f"{c.n_boundary_values} values equal to 0 or 1.", ""]
        for _, r in mm.iterrows():
            lines.append(f"- {r.model} [{r.scale}]: {fmt(r.estimate)} (95% CI {fmt(r.ci_low)} to {fmt(r.ci_high)}; {r.ci_method}); p = {fmt_p(r.p)}")
        ws = wsum[wsum.metric == metric].iloc[0]
        lines += ["", f"- Within-benchmark contrast: mean real - simulated = {fmt(ws.mean_diff)} "
                      f"(median {fmt(ws.median_diff)}), negative in {ws.n_negative}/{ws.k_benchmarks} benchmarks; "
                      f"Wilcoxon p = {fmt_p(ws.wilcoxon_p)}, sign test p = {fmt_p(ws.sign_test_p)}.",
                  f"- Same tool within benchmark: mean paired difference {fmt(ws.mean_paired_tool_diff)}, "
                  f"negative in {ws.n_negative_paired}/{ws.k_paired} benchmarks.", ""]
        lb = lobo[(lobo.metric == metric) & (lobo.excluded != "Edgren-derived observations")]
        ed = lobo[(lobo.metric == metric) & (lobo.excluded == "Edgren-derived observations")].iloc[0]
        lines += [f"- Leave-one-benchmark-out (M1, metric scale): {fmt(lb.m1_diff.min())} to {fmt(lb.m1_diff.max())}; "
                  f"(M2): {fmt(lb.m2_diff.min())} to {fmt(lb.m2_diff.max())}.",
                  f"- Excluding Edgren-derived observations: M1 difference {fmt(ed.m1_diff)} (logit p = {fmt_p(ed.m1_logit_p)}); "
                  f"M2 {fmt(ed.m2_diff)} (95% CI {fmt(ed.m2_ci_low)} to {fmt(ed.m2_ci_high)}; p = {fmt_p(ed.m2_p)}).", ""]
    lines += ["All estimates are associations between benchmark settings and reported performance, not causal effects."]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")

    # LaTeX table for the Supplement
    rows = []
    for metric in METRICS:
        mm = models[models.metric == metric]
        boot_rows = mm[mm.ci_method.str.startswith("benchmark")]
        m0 = mm[mm.model.str.startswith("M0")].iloc[0]
        ws = wsum[wsum.metric == metric].iloc[0]
        for _, r in boot_rows.iterrows():
            rows.append([metric, r.model.split(" (")[0], r.scale, fmt(r.estimate), f"{fmt(r.ci_low)} to {fmt(r.ci_high)}",
                         f"{r.n_obs}/{r.n_tools}/{r.n_bio_datasets}/{r.n_benchmarks}"])
        rows.append([metric, "M0 submitted model", "metric", fmt(m0.estimate), f"{fmt(m0.ci_low)} to {fmt(m0.ci_high)} (Wald)",
                     f"{m0.n_obs}/{m0.n_tools}/{m0.n_bio_datasets}/{m0.n_benchmarks}"])
        rows.append([metric, "Within-benchmark contrast", "metric", fmt(ws.mean_diff), f"{ws.n_negative}/{ws.k_benchmarks} negative",
                     f"--/--/--/{ws.k_benchmarks}"])
    table = pd.DataFrame(rows, columns=["Metric", "Model", "Scale", "Real vs simulated", "95% CI", "Obs/tools/datasets/benchmarks"])
    write_latex_table(
        table, REVISION_DIR / "tables" / "supp_crossed_models.tex",
        caption=r"\textbf{Real versus simulated benchmark settings under crossed dependence.} M1: cross-classified "
                r"logit-normal mixed model with random intercepts for tool, biological dataset and benchmark (primary); "
                r"M2: the same model on the metric scale; M0: submitted model, shown for comparison. Confidence "
                r"intervals for M1 and M2 from a cluster bootstrap resampling whole benchmarks. Estimates are associations, not causal effects.",
        label="tab:supp_crossed_models", colspec="lllrll", resize=True)
    print((OUT / "summary.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
