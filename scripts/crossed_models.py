"""Real versus simulated performance under crossed dependence (tool, biological dataset group, benchmark).

Each observation is a tool score in one dataset of one benchmark publication. Tools recur across
benchmarks and the same biological data (above all Edgren) recur across publications, so tool,
biological dataset group (dataset_metadata.biological_dataset_group) and benchmark are crossed.

Real versus simulated (F1, precision, recall; real and simulated observations, spike-in excluded):
  M1 primary      logit(y*) ~ real + (1|tool) + (1|biological dataset group) + (1|benchmark)
                  y* = (y (n - 1) + 0.5) / n (Smithson & Verkuilen 2006), so 0 and 1 are allowed.
                  Reported on the logit scale and as a difference on the metric scale for a typical
                  tool, dataset group and benchmark: expit(b0 + b1) - expit(b0).
  M2 sensitivity  y ~ real + the same random effects (metric scale).
  M3 beta         beta regression with logit link, y* ~ real + tool (fixed). A beta model with
                  crossed random effects (glmmTMB) is not available in Python; dependence is handled
                  by the benchmark cluster bootstrap. Reported as the average difference in predicted means.
  M0 simple       y ~ real + C(benchmark) + (1|tool): benchmark fixed, tool random; treats observations
                  from the same biological dataset as independent. For comparison only.
Confidence intervals: Wald (model-based) and a cluster bootstrap resampling whole benchmarks.
Sensitivity: leave-one-benchmark-out, excluding Edgren-derived observations, within-benchmark contrasts.

Real observations are also split into cell lines, patient-derived material and mixed sets
(dataset_metadata.real_subtype), each compared with simulated data in the M1/M2 models.
New-tool versus independent benchmarks is reported descriptively, without p-values.

All estimates are associations between benchmark settings and reported performance, not causal effects.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / ".mplconfig"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.othermod.betareg import BetaModel

from new_tool_benchmark_performance import canonical_family, parse_new_tool
from real_simulated_sensitivity import load_table, normalize_columns

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "raw" / "master_analysis_input.xlsx"
OUT = ROOT / "results" / "crossed_models"
METRICS = ["F1", "Precision", "Recall"]
VC = {"tool": "0 + C(Tool)", "bio": "0 + C(Bio)", "bench": "0 + C(Benchmark)"}
SUBTYPES = ["cell line", "patient-derived", "mixed"]


def expit(x):
    return 1.0 / (1.0 + np.exp(-x))


def load() -> pd.DataFrame:
    df = normalize_columns(load_table(INPUT))
    dm = pd.read_excel(INPUT, sheet_name="dataset_metadata")[["dataset_id", "biological_dataset_group", "real_subtype"]]
    df = df.merge(dm, left_on="Dataset", right_on="dataset_id", how="left", validate="many_to_one")
    assert df["biological_dataset_group"].notna().all()
    df["Bio"] = df["biological_dataset_group"].astype(str)
    df["Benchmark"] = df["Benchmark"].astype(str)
    df["real"] = df["DatasetType"].eq("real").astype(int)
    df["Subtype"] = np.where(df["real"].eq(1), df["real_subtype"], "simulated")
    assert df.loc[df["real"].eq(1), "Subtype"].isin(SUBTYPES).all()
    return df


def prepare(df: pd.DataFrame, metric: str) -> pd.DataFrame:
    d = df.dropna(subset=[metric]).copy()
    d["y"] = d[metric].astype(float).clip(0, 1)
    n = len(d)
    d["ys"] = (d["y"] * (n - 1) + 0.5) / n
    d["yl"] = np.log(d["ys"] / (1 - d["ys"]))
    return d


def counts(d: pd.DataFrame) -> dict[str, int]:
    return {"n_obs": len(d), "n_tools": d["Tool"].nunique(), "n_bio_groups": d["Bio"].nunique(),
            "n_datasets": d["Dataset"].nunique(), "n_benchmarks": d["Benchmark"].nunique(),
            "n_real": int(d["real"].sum()), "n_simulated": int((1 - d["real"]).sum()),
            "n_boundary_0_or_1": int(((d["y"] == 0) | (d["y"] == 1)).sum())}


def fit_crossed(d: pd.DataFrame, formula: str):
    """Crossed random-intercept model (REML).

    If the default optimizer does not converge, the model is refitted with other optimizers and
    the converged fit with the highest likelihood is kept; `converged` of the result says whether
    any optimizer converged.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        model = smf.mixedlm(formula, d, groups=np.ones(len(d)), re_formula="0", vc_formula=VC)
        best = model.fit(reml=True)
        if best.converged:
            return best
        for method in (["lbfgs"], ["powell"], ["nm", "lbfgs"], ["cg"]):
            try:
                alt = model.fit(reml=True, method=method)
            except Exception:
                continue
            if alt.converged and (not best.converged or alt.llf > best.llf):
                best = alt
        return best


def fit_beta(d: pd.DataFrame):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return BetaModel.from_formula("ys ~ real + C(Tool)", d).fit(disp=0)


def beta_diff(m3) -> float:
    k = m3.model.exog.shape[1]
    i = list(m3.params.index[:k]).index("real")
    x1, x0 = m3.model.exog.copy(), m3.model.exog.copy()
    x1[:, i], x0[:, i] = 1, 0
    beta = np.asarray(m3.params)[:k]
    return float(np.mean(expit(x1 @ beta) - expit(x0 @ beta)))


def estimates(d: pd.DataFrame, beta: bool = True) -> dict[str, float]:
    m1, m2 = fit_crossed(d, "yl ~ real"), fit_crossed(d, "y ~ real")
    b0, b1 = m1.params["Intercept"], m1.params["real"]
    out = {"m1_logit": float(b1), "m1_diff": float(expit(b0 + b1) - expit(b0)), "m2_diff": float(m2.params["real"]),
           "m1_converged": bool(m1.converged), "m2_converged": bool(m2.converged)}
    if beta:
        m3 = fit_beta(d)
        out.update({"m3_logit": float(m3.params["real"]), "m3_diff": beta_diff(m3)})
    return out


_BOOT_DATA: pd.DataFrame | None = None


def _init_worker(d: pd.DataFrame) -> None:
    global _BOOT_DATA
    _BOOT_DATA = d


def _boot(seed: int):
    d = _BOOT_DATA
    rng = np.random.default_rng(seed)
    benches = d["Benchmark"].unique()
    pick = rng.choice(benches, size=len(benches), replace=True)
    s = pd.concat([d[d["Benchmark"] == b].assign(Benchmark=f"{b}#{i}") for i, b in enumerate(pick)], ignore_index=True)
    if s["real"].nunique() < 2:
        return None
    try:
        return estimates(s)
    except Exception:
        return None


def cluster_bootstrap(d: pd.DataFrame, n_boot: int, seed: int, workers: int) -> pd.DataFrame:
    seeds = [seed + i for i in range(n_boot)]
    if workers > 1:
        # spawn, not fork: forked workers deadlock with the macOS numerical libraries
        with mp.get_context("spawn").Pool(workers, initializer=_init_worker, initargs=(d,)) as pool:
            res = pool.map(_boot, seeds, chunksize=8)
    else:
        _init_worker(d)
        res = [_boot(s) for s in seeds]
    return pd.DataFrame([r for r in res if r is not None])


def within_benchmark(d: pd.DataFrame, metric: str) -> pd.DataFrame:
    rows = []
    for bench, g in d.groupby("Benchmark"):
        if g["real"].nunique() < 2:
            continue
        both = sorted(set(g.loc[g.real == 0, "Tool"]) & set(g.loc[g.real == 1, "Tool"]))
        paired = [g.loc[(g.Tool == t) & (g.real == 1), "y"].mean() - g.loc[(g.Tool == t) & (g.real == 0), "y"].mean() for t in both]
        rows.append({"metric": metric, "benchmark": bench, "n_sim": int((g.real == 0).sum()), "n_real": int((g.real == 1).sum()),
                     "mean_sim": g.loc[g.real == 0, "y"].mean(), "mean_real": g.loc[g.real == 1, "y"].mean(),
                     "diff_real_minus_sim": g.loc[g.real == 1, "y"].mean() - g.loc[g.real == 0, "y"].mean(),
                     "n_tools_both": len(both), "paired_tool_diff_mean": float(np.mean(paired)) if paired else np.nan})
    return pd.DataFrame(rows)


def subtype_models(d: pd.DataFrame, metric: str) -> list[dict]:
    rows = []
    formula = "{} ~ C(Subtype, Treatment('simulated'))"
    m1, m2 = fit_crossed(d, formula.format("yl")), fit_crossed(d, formula.format("y"))
    b0 = m1.params["Intercept"]
    for level in SUBTYPES:
        term = f"C(Subtype, Treatment('simulated'))[T.{level}]"
        if term not in m1.params:
            continue
        sub = d[d.Subtype == level]
        lo, hi = m1.conf_int().loc[term]
        rows.append({"metric": metric, "comparison": f"{level} vs simulated",
                     "n_obs": len(sub), "n_datasets": sub.Dataset.nunique(), "n_bio_groups": sub.Bio.nunique(),
                     "n_benchmarks": sub.Benchmark.nunique(), "mean": sub.y.mean(), "median": sub.y.median(),
                     "m1_logit": m1.params[term], "m1_logit_ci_low": lo, "m1_logit_ci_high": hi, "m1_p": m1.pvalues[term],
                     "m1_diff": expit(b0 + m1.params[term]) - expit(b0),
                     "m1_diff_ci_low": expit(b0 + lo) - expit(b0), "m1_diff_ci_high": expit(b0 + hi) - expit(b0),
                     "m2_diff": m2.params[term], "m2_ci_low": m2.conf_int().loc[term][0], "m2_ci_high": m2.conf_int().loc[term][1],
                     "m2_p": m2.pvalues[term]})
    tp = "C(Subtype, Treatment('simulated'))[T.patient-derived]"
    tc = "C(Subtype, Treatment('simulated'))[T.cell line]"
    if tp in m1.params and tc in m1.params:
        def contrast(m):
            names = list(m.model.exog_names)  # fixed effects only
            r = np.zeros((1, len(names)))
            r[0, names.index(tp)], r[0, names.index(tc)] = 1, -1
            return m.t_test(r)

        t1, t2 = contrast(m1), contrast(m2)
        ci1, ci2 = np.asarray(t1.conf_int()).ravel(), np.asarray(t2.conf_int()).ravel()
        rows.append({"metric": metric, "comparison": "patient-derived vs cell line",
                     "m1_logit": float(np.asarray(t1.effect).ravel()[0]), "m1_logit_ci_low": ci1[0], "m1_logit_ci_high": ci1[1],
                     "m1_p": float(np.asarray(t1.pvalue).ravel()[0]),
                     "m2_diff": float(np.asarray(t2.effect).ravel()[0]), "m2_ci_low": ci2[0], "m2_ci_high": ci2[1],
                     "m2_p": float(np.asarray(t2.pvalue).ravel()[0])})
    return rows


def purpose_descriptive(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = df.dropna(subset=["Benchmark_purpose"]).copy()
    own = d["new_tool"].map(lambda v: {canonical_family(t) for t in parse_new_tool(v)})
    d["is_own_new_tool"] = [canonical_family(t) in fams for t, fams in zip(d["Tool"], own)]
    rows, tool_rows = [], []
    for metric in METRICS:
        for subset, s in [("all tools", d), ("without the authors' own new tool", d[~d["is_own_new_tool"]])]:
            s = s.dropna(subset=[metric])
            for purpose, g in s.groupby("Benchmark_purpose"):
                q1, q3 = g[metric].quantile([0.25, 0.75])
                rows.append({"metric": metric, "subset": subset, "purpose": purpose, "n_obs": len(g),
                             "n_benchmarks": g["Benchmark"].nunique(), "n_tools": g["Tool"].nunique(),
                             "mean": g[metric].mean(), "median": g[metric].median(), "q1": q1, "q3": q3})
            means = s.groupby(["Tool", "Benchmark_purpose"])[metric].mean().unstack()
            if {"new_tool", "independent_comparison"} <= set(means.columns):
                both = means.dropna(subset=["new_tool", "independent_comparison"])
                diff = both["new_tool"] - both["independent_comparison"]
                tool_rows.append({"metric": metric, "subset": subset, "n_tools_in_both": len(both),
                                  "n_higher_in_new_tool": int((diff > 0).sum()),
                                  "median_tool_diff_new_minus_independent": float(diff.median()) if len(diff) else np.nan,
                                  "mean_tool_diff_new_minus_independent": float(diff.mean()) if len(diff) else np.nan})
    return pd.DataFrame(rows), pd.DataFrame(tool_rows)


def fmt(x, d=3):
    return "–" if pd.isna(x) else f"{x:.{d}f}"


def fmt_p(p):
    return "–" if pd.isna(p) else (f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}")


def forest(models: pd.DataFrame, subtypes: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    labels = {"M0": "M0 benchmark fixed, tool random (Wald)", "M1": "M1 crossed logit (bootstrap)", "M2": "M2 crossed linear (bootstrap)",
              "M3": "M3 beta regression (bootstrap)"}
    colors = {"M0": "#9e9e9e", "M1": "#1f4e79", "M2": "#4f81bd", "M3": "#c0504d"}
    ax = axes[0]
    y, ticks, ticklabels = 0, [], []
    for metric in METRICS:
        for key in ["M0", "M1", "M2", "M3"]:
            r = models[(models.metric == metric) & (models.model == key) & (models.scale == "metric difference")]
            r = r[r.ci_method.str.startswith("Wald")] if key == "M0" else r[r.ci_method.str.startswith("benchmark")]
            if r.empty:
                continue
            r = r.iloc[0]
            ax.errorbar(r.estimate, y, xerr=[[r.estimate - r.ci_low], [r.ci_high - r.estimate]], fmt="o", color=colors[key],
                        label=labels[key] if metric == METRICS[0] else None, capsize=3)
            y -= 1
        ticks.append(y + 2.5)
        ticklabels.append(metric)
        y -= 1
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(ticks, ticklabels)
    ax.set_xlabel("Real − simulated (metric scale)")
    ax.set_title("Real versus simulated")
    ax.legend(fontsize=8, loc="lower left")
    ax = axes[1]
    y, ticks, ticklabels = 0, [], []
    sub_colors = {"cell line vs simulated": "#1f4e79", "patient-derived vs simulated": "#c0504d", "mixed vs simulated": "#9bbb59"}
    for metric in METRICS:
        for comp, color in sub_colors.items():
            r = subtypes[(subtypes.metric == metric) & (subtypes.comparison == comp)]
            if r.empty:
                continue
            r = r.iloc[0]
            ax.errorbar(r.m1_diff, y, xerr=[[r.m1_diff - r.m1_diff_ci_low], [r.m1_diff_ci_high - r.m1_diff]], fmt="o",
                        color=color, capsize=3, label=f"{comp.split(' vs')[0]} (n={int(r.n_obs)})" if metric == METRICS[0] else None)
            y -= 1
        ticks.append(y + 2)
        ticklabels.append(metric)
        y -= 1
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(ticks, ticklabels)
    ax.set_xlabel("Real subtype − simulated (M1, metric scale, Wald CI)")
    ax.set_title("Real data by material")
    ax.legend(fontsize=8, loc="lower left")
    fig.suptitle("Crossed models: tool × biological dataset group × benchmark", fontweight="bold")
    fig.tight_layout()
    fig.savefig(path.with_suffix(".png"), dpi=300)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--n-boot", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20261122)
    parser.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    df = load()
    model_rows, sens_rows, within_rows, within_sum, vc_rows, sub_rows = [], [], [], [], [], []
    boot_converged: dict[str, tuple[int, int, int]] = {}
    for metric in METRICS:
        d = prepare(df, metric)
        c = counts(d)
        k_both = int(d.groupby("Benchmark")["real"].nunique().eq(2).sum())
        base = {"metric": metric, **c, "k_benchmarks_both_types": k_both}

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m0 = smf.mixedlm("y ~ real + C(Benchmark)", d, groups=d["Tool"]).fit(reml=True)
        lo, hi = m0.conf_int().loc["real"]
        model_rows.append({**base, "model": "M0", "description": "benchmark fixed, tool random (ignores dataset dependence)", "scale": "metric difference",
                           "estimate": m0.params["real"], "ci_low": lo, "ci_high": hi, "ci_method": "Wald", "p": m0.pvalues["real"]})

        m1, m2, m3 = fit_crossed(d, "yl ~ real"), fit_crossed(d, "y ~ real"), fit_beta(d)
        est = estimates(d)
        boot = cluster_bootstrap(d, args.n_boot, args.seed, args.workers)
        boot.to_csv(OUT / f"bootstrap_{metric.lower()}.csv", index=False)
        boot_converged[metric] = (int(boot["m1_converged"].sum()), int(boot["m2_converged"].sum()), len(boot))

        def bq(col):
            return np.nanpercentile(boot[col], [2.5, 97.5]) if len(boot) else (np.nan, np.nan)

        def bp(col):
            return float(2 * min((boot[col] >= 0).mean(), (boot[col] <= 0).mean())) if len(boot) else np.nan

        b0 = m1.params["Intercept"]
        w1 = m1.conf_int().loc["real"]
        spec = [
            ("M1", "crossed logit-normal (primary)", "logit", est["m1_logit"], w1, m1.pvalues["real"], "m1_logit"),
            ("M1", "crossed logit-normal (primary)", "metric difference", est["m1_diff"],
             (expit(b0 + w1[0]) - expit(b0), expit(b0 + w1[1]) - expit(b0)), m1.pvalues["real"], "m1_diff"),
            ("M2", "crossed linear (sensitivity)", "metric difference", est["m2_diff"], m2.conf_int().loc["real"], m2.pvalues["real"], "m2_diff"),
            ("M3", "beta regression, tool fixed", "logit", est["m3_logit"], m3.conf_int().loc["real"], m3.pvalues["real"], "m3_logit"),
            ("M3", "beta regression, tool fixed", "metric difference", est["m3_diff"], (np.nan, np.nan), np.nan, "m3_diff"),
        ]
        for key, desc, scale, value, wald, p, col in spec:
            row = {**base, "model": key, "description": desc, "scale": scale, "estimate": value}
            lo, hi = bq(col)
            model_rows.append({**row, "ci_low": lo, "ci_high": hi, "ci_method": f"benchmark cluster bootstrap (B={len(boot)})",
                               "p": bp(col)})
            if np.isfinite(wald[0]):
                model_rows.append({**row, "ci_low": wald[0], "ci_high": wald[1], "ci_method": "Wald", "p": p})

        vc_rows.append({"metric": metric, "model": "M1 (logit)", **dict(zip(m1.model.exog_vc.names, m1.vcomp)), "residual": m1.scale})
        vc_rows.append({"metric": metric, "model": "M2 (metric)", **dict(zip(m2.model.exog_vc.names, m2.vcomp)), "residual": m2.scale})

        for bench in sorted(d["Benchmark"].unique(), key=int):
            s = d[d["Benchmark"] != bench]
            if s["real"].nunique() == 2:
                sens_rows.append({"metric": metric, "excluded": f"benchmark {bench}", **estimates(s, beta=False), **counts(s)})
        s = d[d["Is_Edgren_Binary"] != "yes"]
        m1e, m2e = fit_crossed(s, "yl ~ real"), fit_crossed(s, "y ~ real")
        sens_rows.append({"metric": metric, "excluded": "Edgren-derived observations", **estimates(s, beta=False), **counts(s),
                          "m1_logit_p": m1e.pvalues["real"], "m2_p": m2e.pvalues["real"],
                          "m2_ci_low": m2e.conf_int().loc["real"][0], "m2_ci_high": m2e.conf_int().loc["real"][1]})

        w = within_benchmark(d, metric)
        within_rows.append(w)
        diffs = w["diff_real_minus_sim"].to_numpy()
        paired = w["paired_tool_diff_mean"].dropna().to_numpy()
        within_sum.append({"metric": metric, "k_benchmarks": len(diffs), "mean_diff": diffs.mean(), "median_diff": np.median(diffs),
                           "n_negative": int((diffs < 0).sum()),
                           "wilcoxon_p": stats.wilcoxon(diffs).pvalue if len(diffs) >= 5 else np.nan,
                           "sign_test_p": stats.binomtest(int((diffs < 0).sum()), len(diffs), 0.5).pvalue,
                           "k_paired": len(paired), "mean_paired_tool_diff": paired.mean() if len(paired) else np.nan,
                           "n_negative_paired": int((paired < 0).sum())})
        sub_rows += subtype_models(d, metric)

    models, sens = pd.DataFrame(model_rows), pd.DataFrame(sens_rows)
    within, wsum, subtypes = pd.concat(within_rows, ignore_index=True), pd.DataFrame(within_sum), pd.DataFrame(sub_rows)
    purpose, purpose_tools = purpose_descriptive(df)
    models.to_csv(OUT / "real_vs_simulated_models.csv", index=False)
    sens.to_csv(OUT / "sensitivity_lobo_edgren.csv", index=False)
    within.to_csv(OUT / "within_benchmark_contrasts.csv", index=False)
    wsum.to_csv(OUT / "within_benchmark_summary.csv", index=False)
    pd.DataFrame(vc_rows).to_csv(OUT / "variance_components.csv", index=False)
    subtypes.to_csv(OUT / "real_subtype_models.csv", index=False)
    purpose.to_csv(OUT / "benchmark_purpose_descriptive.csv", index=False)
    purpose_tools.to_csv(OUT / "benchmark_purpose_same_tool.csv", index=False)
    forest(models, subtypes, OUT / "crossed_models_forest")

    lines = ["# Crossed-dependence models, beta regression, real-data subtypes", ""]
    for metric in METRICS:
        mm = models[models.metric == metric]
        c = mm.iloc[0]
        lines += [f"## {metric}", "",
                  f"{c.n_obs} observations ({c.n_real} real, {c.n_simulated} simulated), {c.n_tools} tools, "
                  f"{c.n_bio_groups} biological dataset groups, {c.n_benchmarks} benchmarks "
                  f"({c.k_benchmarks_both_types} with both real and simulated data).", ""]
        for _, r in mm.iterrows():
            lines.append(f"- {r.model} {r.description} [{r.scale}]: {fmt(r.estimate)} "
                         f"(95% CI {fmt(r.ci_low)} to {fmt(r.ci_high)}; {r.ci_method}); p = {fmt_p(r.p)}")
        ws = wsum[wsum.metric == metric].iloc[0]
        lb = sens[(sens.metric == metric) & sens.excluded.str.startswith("benchmark")]
        ed = sens[(sens.metric == metric) & (sens.excluded == "Edgren-derived observations")].iloc[0]
        c1, c2, nb = boot_converged[metric]
        lines += ["", f"- Convergence (after optimizer fallback): bootstrap M1 {c1}/{nb}, M2 {c2}/{nb} replicates; "
                      f"leave-one-benchmark-out M1 {int(lb.m1_converged.sum())}/{len(lb)} fits.",
                  f"- Leave-one-benchmark-out, M1 metric difference: {fmt(lb.m1_diff.min())} to {fmt(lb.m1_diff.max())}; "
                      f"M2: {fmt(lb.m2_diff.min())} to {fmt(lb.m2_diff.max())}.",
                  f"- Without Edgren-derived observations: M1 {fmt(ed.m1_diff)} (logit p = {fmt_p(ed.m1_logit_p)}); "
                  f"M2 {fmt(ed.m2_diff)} (95% CI {fmt(ed.m2_ci_low)} to {fmt(ed.m2_ci_high)}; p = {fmt_p(ed.m2_p)}).",
                  f"- Within-benchmark contrast (real − simulated, one value per benchmark): mean {fmt(ws.mean_diff)}, "
                  f"negative in {ws.n_negative}/{ws.k_benchmarks} benchmarks (sign test p = {fmt_p(ws.sign_test_p)}); "
                  f"same tool within benchmark: mean {fmt(ws.mean_paired_tool_diff)}, negative in {ws.n_negative_paired}/{ws.k_paired}.", ""]
        for _, r in subtypes[subtypes.metric == metric].iterrows():
            if r.comparison.endswith("simulated"):
                lines.append(f"- {r.comparison} (n = {int(r.n_obs)}, {int(r.n_datasets)} datasets, {int(r.n_benchmarks)} benchmarks): "
                             f"M1 {fmt(r.m1_diff)} (95% CI {fmt(r.m1_diff_ci_low)} to {fmt(r.m1_diff_ci_high)}; p = {fmt_p(r.m1_p)}); "
                             f"M2 {fmt(r.m2_diff)} (p = {fmt_p(r.m2_p)})")
            else:
                lines.append(f"- {r.comparison}: M1 logit {fmt(r.m1_logit)} (p = {fmt_p(r.m1_p)}); M2 {fmt(r.m2_diff)} "
                             f"(95% CI {fmt(r.m2_ci_low)} to {fmt(r.m2_ci_high)}; p = {fmt_p(r.m2_p)})")
        lines.append("")
    lines += ["## Benchmark purpose (descriptive)", ""]
    for _, r in purpose[purpose.subset == "all tools"].iterrows():
        lines.append(f"- {r.metric}, {r.purpose}: {r.n_obs} observations, {r.n_benchmarks} benchmarks, median {fmt(r["median"])} "
                     f"(IQR {fmt(r.q1)}–{fmt(r.q3)}), mean {fmt(r['mean'])}")
    for _, r in purpose_tools.iterrows():
        lines.append(f"- Same tool, {r.metric} ({r.subset}): higher in new-tool benchmarks for {r.n_higher_in_new_tool}/{r.n_tools_in_both} "
                     f"tools; median difference {fmt(r.median_tool_diff_new_minus_independent)}")
    lines += ["", "All estimates are associations between benchmark settings and reported performance, not causal effects."]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print((OUT / "summary.md").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
