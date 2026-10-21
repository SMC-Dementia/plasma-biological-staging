#!/usr/bin/env python
"""Longitudinal CDR-SB validation (Fig. 3a; Supplementary Table S3)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL, STAGE_LABELS
from biological_staging.io import read_table, require_columns
from biological_staging.models import fit_mixedlm, likelihood_ratio_test, fixed_effect_params_and_cov
from biological_staging.stats import linear_contrast, bh_fdr
from biological_staging.staging import coerce_stage_value

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("--outdir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--id-col",default=COL["id"]); ap.add_argument("--time-col",default=COL["time"]); ap.add_argument("--cdrsb-col",default=COL["cdrsb"]); ap.add_argument("--stage-col",default="stage"); args=ap.parse_args()
    df=read_table(args.input); require_columns(df,[args.id_col,args.time_col,args.cdrsb_col,args.stage_col])
    df=df.copy(); df[args.stage_col]=df[args.stage_col].map(coerce_stage_value); df=df[df[args.stage_col].isin(range(5))].copy(); df[args.stage_col]=df[args.stage_col].astype(int)
    formula=f"Q('{args.cdrsb_col}') ~ {args.time_col} * C({args.stage_col})"
    reduced=f"Q('{args.cdrsb_col}') ~ {args.time_col} + C({args.stage_col})"
    full=fit_mixedlm(formula,df,args.id_col,re_formula=f"~{args.time_col}",reml=False)
    red=fit_mixedlm(reduced,df,args.id_col,re_formula=f"~{args.time_col}",reml=False)
    lr,dfd,p=likelihood_ratio_test(full,red,4)
    params,cov=fixed_effect_params_and_cov(full)

    slopes=[]
    for s in range(5):
        slope_weights={args.time_col:1.0}
        if s>0: slope_weights[f"{args.time_col}:C({args.stage_col})[T.{s}]"]=1.0
        slope=linear_contrast(params,cov,slope_weights)

        intercept_weights={"Intercept":1.0}
        if s>0: intercept_weights[f"C({args.stage_col})[T.{s}]"]=1.0
        intercept=linear_contrast(params,cov,intercept_weights)

        slope.update({
            "stage":s,
            "stage_label":STAGE_LABELS[s],
            "intercept_at_time0":intercept["estimate"],
            "intercept_se":intercept["se"],
            "intercept_ci_low":intercept["ci_low"],
            "intercept_ci_high":intercept["ci_high"],
            "N_participants":df.loc[df[args.stage_col]==s,args.id_col].nunique(),
            "N_assessments":int((df[args.stage_col]==s).sum()),
        })
        slopes.append(slope)
    slopes_df=pd.DataFrame(slopes).rename(columns={"estimate":"annual_slope"})

    adj=[]
    for a,b in zip(range(4),range(1,5)):
        weights={}
        # Difference slope(b)-slope(a); common time term cancels.
        if b>0: weights[f"{args.time_col}:C({args.stage_col})[T.{b}]"]=1.0
        if a>0: weights[f"{args.time_col}:C({args.stage_col})[T.{a}]"]=-1.0
        c=linear_contrast(params,cov,weights); c.update({"stage_from":a,"stage_to":b,"contrast":f"{STAGE_LABELS[a]} → {STAGE_LABELS[b]}"}); adj.append(c)
    adj_df=pd.DataFrame(adj); adj_df["q"]=bh_fdr(adj_df.p)
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)
    slopes_df.to_csv(outdir/"Supplementary_Table_S3A_stage_slopes.csv",index=False)
    adj_df.to_csv(outdir/"Supplementary_Table_S3B_adjacent_slope_contrasts.csv",index=False)
    pd.DataFrame([{"chi2":lr,"df":dfd,"p":p,"participants":df[args.id_col].nunique(),"assessments":len(df)}]).to_csv(outdir/"Figure3A_global_stage_by_time.csv",index=False)
    print(f"Stage×time LR χ²({dfd})={lr:.3f}, P={p:.3g}")

if __name__=="__main__": main()
