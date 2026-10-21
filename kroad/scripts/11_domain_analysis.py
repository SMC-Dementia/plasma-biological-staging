#!/usr/bin/env python
"""Prespecified stage-oriented longitudinal domain scores (Supplementary Table S5C)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL
from biological_staging.io import read_table, canonicalize_columns, require_columns
from biological_staging.models import fit_mixedlm, fixed_effect_params_and_cov
from biological_staging.stats import linear_contrast, bh_fdr
from biological_staging.staging import coerce_stage_value

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("serial_input"); ap.add_argument("baseline_input"); ap.add_argument("--outdir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--id-col",default=COL["id"]); ap.add_argument("--time-col",default=COL["time"]); args=ap.parse_args()
    serial=canonicalize_columns(read_table(args.serial_input)); base=canonicalize_columns(read_table(args.baseline_input)); require_columns(serial,[args.id_col,args.time_col]);
    if "stage" in base: base.stage=base.stage.map(coerce_stage_value); base=base[base.stage.isin(range(5))].copy()
    meta=read_table(ROOT/"metadata"/"stage_associated_34.csv");
    long=serial[[args.id_col,args.time_col]].copy()
    usable=[]
    for _,r in meta.iterrows():
        p=r.protein
        if p not in serial or p not in base: continue
        b=pd.to_numeric(base[p],errors="coerce"); mu=b.mean(); sd=b.std(ddof=1); sign=1 if r.ordinal_beta_reference>=0 else -1
        long[p]=sign*(pd.to_numeric(serial[p],errors="coerce")-mu)/sd; usable.append(p)
    rows=[]
    for domain,g in meta[meta.protein.isin(usable)].groupby("domain",observed=True):
        proteins=g.protein.tolist()
        if len(proteins)<3: continue
        d=long[[args.id_col,args.time_col]+proteins].copy(); d["domain_score"]=d[proteins].mean(axis=1,skipna=True); d=d.dropna(subset=[args.id_col,args.time_col,"domain_score"])
        fit=fit_mixedlm("domain_score ~ "+args.time_col,d,args.id_col,re_formula=f"~{args.time_col}",reml=False); params,cov=fixed_effect_params_and_cov(fit); c=linear_contrast(params,cov,{args.time_col:1.0}); c.update({"domain":domain,"N_proteins":len(proteins),"proteins":"; ".join(proteins),"N_participants":d[args.id_col].nunique(),"N_draws":len(d)}); rows.append(c)
    res=pd.DataFrame(rows).rename(columns={"estimate":"annual_slope_SD_per_year"}); res["q"]=bh_fdr(res.p); outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True); res.to_csv(outdir/"Supplementary_Table_S5C_domain_scores.csv",index=False); print(res[["domain","N_proteins","annual_slope_SD_per_year","q"]].to_string(index=False))
if __name__=="__main__": main()
