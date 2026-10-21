"""Statistical helpers."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Iterable

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve
from statsmodels.stats.multitest import multipletests


def bh_fdr(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    good = np.isfinite(p)
    if good.any():
        q[good] = multipletests(p[good], method="fdr_bh")[1]
    return q


def zscore(series: pd.Series, mean: float | None = None, sd: float | None = None) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    m = x.mean() if mean is None else mean
    s = x.std(ddof=1) if sd is None else sd
    if not np.isfinite(s) or s == 0:
        raise ValueError(f"Cannot standardize {series.name!r}: SD={s}")
    return (x - m) / s


def standardize_within_group(df: pd.DataFrame, cols: list[str], group_col: str) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        out[col] = out.groupby(group_col, observed=True)[col].transform(
            lambda s: (s - s.mean()) / s.std(ddof=1)
        )
    return out


def youden_cutoff(y_true: pd.Series | np.ndarray, score: pd.Series | np.ndarray):
    y = np.asarray(y_true, dtype=int)
    x = np.asarray(score, dtype=float)
    good = np.isfinite(x) & np.isfinite(y)
    fpr, tpr, thresholds = roc_curve(y[good], x[good])
    j = tpr - fpr
    i = int(np.nanargmax(j))
    return {
        "cutoff": float(thresholds[i]),
        "sensitivity": float(tpr[i]),
        "specificity": float(1 - fpr[i]),
        "youden_j": float(j[i]),
    }


def bootstrap_statistic(
    df: pd.DataFrame,
    statistic: Callable[[pd.DataFrame], float],
    *,
    n_boot: int = 5000,
    seed: int = 20260826,
    stratify_col: str | None = None,
):
    """Participant-row bootstrap with optional stratification."""
    rng = np.random.default_rng(seed)
    vals = []
    if stratify_col is None:
        n = len(df)
        for _ in range(n_boot):
            idx = rng.integers(0, n, n)
            try:
                vals.append(float(statistic(df.iloc[idx])))
            except Exception:
                vals.append(np.nan)
    else:
        groups = [g for _, g in df.groupby(stratify_col, observed=True)]
        for _ in range(n_boot):
            sample = []
            for g in groups:
                idx = rng.integers(0, len(g), len(g))
                sample.append(g.iloc[idx])
            b = pd.concat(sample, ignore_index=True)
            try:
                vals.append(float(statistic(b)))
            except Exception:
                vals.append(np.nan)
    a = np.asarray(vals, dtype=float)
    a = a[np.isfinite(a)]
    return a


def auc_ci_bootstrap(y, score, *, n_boot=5000, seed=20260826):
    d = pd.DataFrame({"y": y, "x": score}).dropna()
    point = roc_auc_score(d["y"], d["x"])
    boots = bootstrap_statistic(
        d,
        lambda b: roc_auc_score(b["y"], b["x"]) if b["y"].nunique() == 2 else np.nan,
        n_boot=n_boot,
        seed=seed,
        stratify_col="y",
    )
    return point, tuple(np.quantile(boots, [0.025, 0.975]))


def paired_auc_difference_bootstrap(y, score_a, score_b, *, n_boot=5000, seed=20260826):
    """Bootstrap paired ΔAUC = AUC(score_a) - AUC(score_b)."""
    d = pd.DataFrame({"y": y, "a": score_a, "b": score_b}).dropna()
    point = roc_auc_score(d.y, d.a) - roc_auc_score(d.y, d.b)

    def stat_fn(b):
        if b.y.nunique() < 2:
            return np.nan
        return roc_auc_score(b.y, b.a) - roc_auc_score(b.y, b.b)

    boots = bootstrap_statistic(d, stat_fn, n_boot=n_boot, seed=seed, stratify_col="y")
    ci = tuple(np.quantile(boots, [0.025, 0.975]))
    # Two-sided bootstrap sign test; floor avoids exact zero.
    p = 2 * min((boots <= 0).mean(), (boots >= 0).mean())
    p = min(1.0, float(p))
    return float(point), ci, p


def bootstrap_youden_cutoff(y, score, *, n_boot=5000, seed=20260826):
    d = pd.DataFrame({"y": y, "x": score}).dropna()
    point = youden_cutoff(d.y, d.x)
    boots = bootstrap_statistic(
        d,
        lambda b: youden_cutoff(b.y, b.x)["cutoff"] if b.y.nunique() == 2 else np.nan,
        n_boot=n_boot,
        seed=seed,
        stratify_col="y",
    )
    point["cutoff_ci_low"], point["cutoff_ci_high"] = np.quantile(boots, [0.025, 0.975])
    return point


def spearman_with_bootstrap_ci(x, y, *, n_boot=5000, seed=20260826):
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    rho, p = stats.spearmanr(d.x, d.y)
    boots = bootstrap_statistic(
        d,
        lambda b: stats.spearmanr(b.x, b.y).statistic,
        n_boot=n_boot,
        seed=seed,
    )
    lo, hi = np.quantile(boots, [0.025, 0.975])
    return float(rho), float(p), (float(lo), float(hi))


def partial_r2_from_sse(sse_reduced: float, sse_full: float) -> float:
    """Incremental partial R² = (SSE_reduced - SSE_full) / SSE_reduced."""
    if sse_reduced <= 0:
        return np.nan
    return float((sse_reduced - sse_full) / sse_reduced)


def linear_contrast(params: pd.Series, cov: pd.DataFrame, weights: dict[str, float]):
    """Wald contrast from fitted fixed-effect parameters."""
    names = list(params.index)
    missing=set(weights)-set(names)
    if missing:
        raise KeyError(f'Contrast references absent model terms: {sorted(missing)}')
    w = np.array([weights.get(n, 0.0) for n in names], dtype=float)
    est = float(w @ params.to_numpy())
    v = float(w @ cov.loc[names, names].to_numpy() @ w)
    se = float(np.sqrt(max(v, 0)))
    z = est / se if se > 0 else np.nan
    p = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    return {"estimate": est, "se": se, "z": z, "p": p, "ci_low": est - 1.96 * se, "ci_high": est + 1.96 * se}


def exact_directional_binomial(n_progression: int, n_regression: int):
    n = int(n_progression + n_regression)
    if n == 0:
        return np.nan
    return float(stats.binomtest(n_progression, n=n, p=0.5, alternative="two-sided").pvalue)
