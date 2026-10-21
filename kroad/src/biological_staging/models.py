"""Regression and mixed-model wrappers."""
from __future__ import annotations

import warnings
import os
import json
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

from .stats import linear_contrast, partial_r2_from_sse


def fit_ols(formula: str, data: pd.DataFrame, robust: str | None = None):
    model = smf.ols(formula, data=data).fit()
    return model.get_robustcov_results(cov_type=robust) if robust else model


def fit_mixedlm(
    formula: str,
    data: pd.DataFrame,
    group_col: str,
    re_formula: str = "~years_from_baseline",
    *,
    reml: bool = False,
    method: str = "lbfgs",
    maxiter: int = 2000,
):
    """Random intercept and slope model. Alternative optimizers are tried if the first does not converge."""
    attempts=[]
    candidates=[]
    for optimizer in dict.fromkeys([method, 'bfgs', 'powell']):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            try:
                md = smf.mixedlm(formula, data=data, groups=data[group_col], re_formula=re_formula)
                result = md.fit(reml=reml, method=optimizer, maxiter=maxiter, disp=False)
                ok=bool(result.converged) and np.isfinite(result.llf)
                attempts.append({'optimizer':optimizer,'converged':bool(ok),'llf':float(result.llf),'warnings':[str(x.message) for x in caught]})
                if ok: candidates.append(result)
            except Exception as exc:
                attempts.append({'optimizer':optimizer,'converged':False,'error':str(exc)})
        if candidates: break
    log=os.environ.get('BIOSTAGING_MODEL_LOG')
    if log:
        with open(log,'a',encoding='utf-8') as stream:
            stream.write(json.dumps({'formula':formula,'n':len(data),'groups':data[group_col].nunique(),'random_effects':re_formula,'attempts':attempts},default=str)+'\n')
    if not candidates:
        raise RuntimeError(f'No converged random intercept and slope fit: {formula}; {attempts}')
    fit=candidates[0]
    if len(attempts)>1:
        warnings.warn(f'MixedLM required optimizer retry: {formula}; selected {attempts[-1]["optimizer"]}')
    return fit


def likelihood_ratio_test(full_fit, reduced_fit, df_diff: int):
    stat = 2 * (full_fit.llf - reduced_fit.llf)
    p = stats.chi2.sf(stat, df_diff)
    return float(stat), int(df_diff), float(p)


def categorical_stage_omnibus(data: pd.DataFrame, outcome: str, age: str, sex: str, stage: str = "stage"):
    """Primary cross-sectional model: protein ~ C(stage) + age + sex."""
    full = smf.ols(f"Q('{outcome}') ~ C({stage}) + {age} + C({sex})", data=data).fit()
    reduced = smf.ols(f"Q('{outcome}') ~ {age} + C({sex})", data=data).fit()
    # Nested-model F test.
    f_stat, p, df_diff = full.compare_f_test(reduced)
    partial_r2 = partial_r2_from_sse(np.sum(reduced.resid**2), np.sum(full.resid**2))
    return full, {
        "stage_f": float(f_stat),
        "omnibus_p": float(p),
        "df_diff": int(df_diff),
        "partial_r2": float(partial_r2),
    }


def ordinal_stage_robust(data: pd.DataFrame, outcome: str, age: str, sex: str, stage: str = "stage"):
    fit = smf.ols(f"Q('{outcome}') ~ {stage} + {age} + C({sex})", data=data).fit(cov_type="HC3")
    return {
        "beta": float(fit.params[stage]),
        "se": float(fit.bse[stage]),
        "p": float(fit.pvalues[stage]),
        "fit": fit,
    }


def adjacent_stage_contrasts_robust(data: pd.DataFrame, outcome: str, age: str, sex: str, stage: str = "stage"):
    fit = smf.ols(f"Q('{outcome}') ~ C({stage}) + {age} + C({sex})", data=data).fit(cov_type="HC3")
    params = fit.params
    cov = fit.cov_params()

    def term(s):
        return None if s == 0 else f"C({stage})[T.{s}]"

    rows = []
    for a, b in zip(range(4), range(1, 5)):
        weights = {}
        ta, tb = term(a), term(b)
        if tb is not None:
            weights[tb] = 1.0
        if ta is not None:
            weights[ta] = -1.0
        c = linear_contrast(params, cov, weights)
        c.update({"stage_from": a, "stage_to": b})
        rows.append(c)
    return fit, pd.DataFrame(rows)


def fixed_effect_params_and_cov(fit):
    """Extract only fixed effects from a statsmodels MixedLM result."""
    names = list(fit.fe_params.index)
    cov = fit.cov_params().loc[names, names]
    return fit.fe_params, cov
