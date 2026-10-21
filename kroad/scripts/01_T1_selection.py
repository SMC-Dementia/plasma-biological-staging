#!/usr/bin/env python
"""T1 marker benchmarking for Aβ PET positivity (Fig. 2a, Suppl. Fig. S1/Table S1)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from biological_staging.config import AMYLOID_CENTILOID_POSITIVE, COL
from biological_staging.io import read_table, canonicalize_columns
from biological_staging.stats import (
    auc_ci_bootstrap, paired_auc_difference_bootstrap, bootstrap_youden_cutoff,
    spearman_with_bootstrap_ci, standardize_within_group
)

PAIRS = [
    ("p-tau181", "BD-p-tau181"),
    ("p-tau212", "BD-p-tau212"),
    ("p-tau217", "BD-p-tau217"),
    ("p-tau231", "BD-p-tau231"),
]
BD_EXPANDED = ["BD-p-tau181", "BD-p-tau212", "BD-p-tau217", "BD-p-tau231"]
BD_CANONICAL = ["BD-p-tau181", "BD-p-tau217", "BD-p-tau231"]
PTAU_EXPANDED = ["p-tau181", "p-tau212", "p-tau217", "p-tau231"]
PTAU_CANONICAL = ["p-tau181", "p-tau217", "p-tau231"]


def equal_weight_composite(df, cols, phase_col=None):
    d = df.copy()
    if phase_col and phase_col in d.columns:
        d = standardize_within_group(d, cols, phase_col)
    else:
        for c in cols:
            d[c] = (d[c] - d[c].mean()) / d[c].std(ddof=1)
    return d[cols].mean(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("--outdir", default=str(ROOT / "outputs" / "tables"))
    ap.add_argument("--centiloid-col", default=COL["centiloid"])
    ap.add_argument("--phase-col", default=COL["phase"])
    ap.add_argument("--n-boot", type=int, default=2000, help="Bootstrap replicates (2,000 in the manuscript)")
    ap.add_argument("--seed", type=int, default=20260826)
    args = ap.parse_args()

    df = canonicalize_columns(read_table(args.input))
    cl = pd.to_numeric(df[args.centiloid_col], errors="coerce")
    df["amyloid_positive"] = (cl > AMYLOID_CENTILOID_POSITIVE).where(cl.notna()).astype(float)
    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    for _, marker in PAIRS:
        if marker not in df: continue
        d = df[["amyloid_positive", args.centiloid_col, marker]].dropna()
        auc, auc_ci = auc_ci_bootstrap(d.amyloid_positive, d[marker], n_boot=args.n_boot, seed=args.seed)
        rho, _, rho_ci = spearman_with_bootstrap_ci(d[marker], d[args.centiloid_col], n_boot=args.n_boot, seed=args.seed)
        rows.append({"assay": marker, "N": len(d), "AUC": auc, "AUC_CI_low": auc_ci[0], "AUC_CI_high": auc_ci[1], "Spearman_rho": rho, "rho_CI_low": rho_ci[0], "rho_CI_high": rho_ci[1]})
        conventional = marker.replace("BD-", "")
        if conventional in df:
            d2 = df[["amyloid_positive", marker, conventional]].dropna()
            auc2, ci2 = auc_ci_bootstrap(d2.amyloid_positive, d2[conventional], n_boot=args.n_boot, seed=args.seed)
            rho2, _, rho_ci2 = spearman_with_bootstrap_ci(d2[conventional], df.loc[d2.index, args.centiloid_col], n_boot=args.n_boot, seed=args.seed)
            rows.append({"assay": conventional, "N": len(d2), "AUC": auc2, "AUC_CI_low": ci2[0], "AUC_CI_high": ci2[1], "Spearman_rho": rho2, "rho_CI_low": rho_ci2[0], "rho_CI_high": rho_ci2[1]})
    pd.DataFrame(rows).drop_duplicates("assay").to_csv(outdir / "Supplementary_Table_S1A_T1_candidates.csv", index=False)

    pair_rows = []
    for conventional, bd in PAIRS:
        if conventional not in df or bd not in df: continue
        d = df[["amyloid_positive", conventional, bd]].dropna()
        delta, ci, p = paired_auc_difference_bootstrap(d.amyloid_positive, d[bd], d[conventional], n_boot=args.n_boot, seed=args.seed)
        pair_rows.append({"phosphosite": conventional, "assay_pair_rho": d[[conventional, bd]].corr(method="spearman").iloc[0,1], "BD_AUC": roc_auc_score(d.amyloid_positive, d[bd]), "p_tau_AUC": roc_auc_score(d.amyloid_positive, d[conventional]), "delta_AUC_BD_minus_p": delta, "CI_low": ci[0], "CI_high": ci[1], "bootstrap_p": p})
    pd.DataFrame(pair_rows).to_csv(outdir / "Supplementary_Table_S1B_paired_AUC.csv", index=False)

    representations = {"BD-p-tau217 alone": df["BD-p-tau217"]}
    for name, cols in {
        "BD canonical (181+217+231)": BD_CANONICAL,
        "BD expanded (+212)": BD_EXPANDED,
        "p-tau canonical (181+217+231)": PTAU_CANONICAL,
        "p-tau expanded (+212)": PTAU_EXPANDED,
    }.items():
        if all(c in df.columns for c in cols):
            representations[name] = equal_weight_composite(df, cols, args.phase_col if args.phase_col in df.columns else None)

    # Site-consensus representations: within-phase standardised plain and BD assays
    # at each phosphosite, averaged with equal weight.
    phase_col = args.phase_col if args.phase_col in df.columns else None
    site_pairs = {
        "181": ["p-tau181", "BD-p-tau181"],
        "212": ["p-tau212", "BD-p-tau212"],
        "217": ["p-tau217", "BD-p-tau217"],
        "231": ["p-tau231", "BD-p-tau231"],
    }
    if all(all(c in df.columns for c in site_pairs[k]) for k in ["181","217","231"]):
        zcols = sum([site_pairs[k] for k in ["181","217","231"]], [])
        representations["Site-consensus canonical"] = equal_weight_composite(df, zcols, phase_col)
    if all(all(c in df.columns for c in site_pairs[k]) for k in ["181","212","217","231"]):
        zcols = sum([site_pairs[k] for k in ["181","212","217","231"]], [])
        representations["Site-consensus expanded (+212)"] = equal_weight_composite(df, zcols, phase_col)
    if all(c in df.columns for c in site_pairs["217"]):
        representations["p-tau217 site-consensus"] = equal_weight_composite(df, site_pairs["217"], phase_col)

    rep_rows = []
    ref = representations["BD-p-tau217 alone"]
    for name, score in representations.items():
        d = pd.DataFrame({"y": df.amyloid_positive, "ref": ref, "score": score, "centiloid": df[args.centiloid_col]}).dropna()
        auc, ci = auc_ci_bootstrap(d.y, d.score, n_boot=args.n_boot, seed=args.seed)
        rho = d[["score", "centiloid"]].corr(method="spearman").iloc[0,1]
        delta, dci, p = paired_auc_difference_bootstrap(d.y, d.ref, d.score, n_boot=args.n_boot, seed=args.seed)
        rep_rows.append({"representation": name, "AUC": auc, "CI_low": ci[0], "CI_high": ci[1], "Spearman_rho": rho, "delta_AUC_BD217_minus_model": delta, "delta_CI_low": dci[0], "delta_CI_high": dci[1], "bootstrap_p": p})
    pd.DataFrame(rep_rows).to_csv(outdir / "Supplementary_Table_S1C_multianalyte.csv", index=False)

    d = df[["amyloid_positive", "BD-p-tau217"]].dropna()
    cutoff = bootstrap_youden_cutoff(d.amyloid_positive, d["BD-p-tau217"], n_boot=args.n_boot, seed=args.seed)
    pd.DataFrame([{"N":len(d), "amyloid_positive_N":int(d.amyloid_positive.sum()), **cutoff, "AUC":roc_auc_score(d.amyloid_positive,d["BD-p-tau217"])}]).to_csv(outdir / "Supplementary_Table_S1D_T1_threshold.csv", index=False)
    print(f"T1 outputs written to {outdir}")

if __name__ == "__main__": main()
