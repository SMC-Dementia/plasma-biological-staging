#!/usr/bin/env python
"""Longitudinal recapitulation of the 34 stage-associated proteins (Fig. 4d; Suppl. S3/S5B)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
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
    if "stage" in base: base["stage"]=base.stage.map(coerce_stage_value); base=base[base.stage.isin(range(5))].copy()
    meta=read_table(ROOT/"metadata"/"stage_associated_34.csv"); proteins=[p for p in meta.protein if p in serial.columns and p in base.columns]
    rows=[]
    for p in proteins:
        bx=pd.to_numeric(base[p],errors="coerce"); mu=float(bx.mean()); sd=float(bx.std(ddof=1));
        d=serial[[args.id_col,args.time_col,p]].copy(); d["protein_z"]=(pd.to_numeric(d[p],errors="coerce")-mu)/sd; d=d.dropna()
        fit=fit_mixedlm("protein_z ~ "+args.time_col,d,args.id_col,re_formula=f"~{args.time_col}",reml=False)
        params,cov=fixed_effect_params_and_cov(fit); c=linear_contrast(params,cov,{args.time_col:1.0})
        beta=float(meta.loc[meta.protein==p,"ordinal_beta_reference"].iloc[0]); orientation="Increasing" if beta>=0 else "Decreasing"; sign=1.0 if beta>=0 else -1.0
        c.update({"protein":p,"domain":meta.loc[meta.protein==p,"domain"].iloc[0],"cross_sectional_stage_orientation":orientation,"orientation_sign":sign,"stage_oriented_slope":sign*c["estimate"],"N_participants":d[args.id_col].nunique(),"N_draws":len(d),"baseline_mean":mu,"baseline_sd":sd}); rows.append(c)
    res=pd.DataFrame(rows).rename(columns={"estimate":"annual_slope_SD_per_year"}); res["q"]=bh_fdr(res.p); res["stage_consistent"]=(res.stage_oriented_slope>0); res["fdr_significant"]=res.q<0.05
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True); res=res.sort_values("q"); res.to_csv(outdir/"Supplementary_Table_S5B_longitudinal_34.csv",index=False); res[res.fdr_significant].to_csv(outdir/"Supplementary_Figure_S3_longitudinal_FDR15.csv",index=False)
    sig=res[res.fdr_significant]; print(f"Longitudinal q<0.05: {len(sig)}; stage-consistent: {int(sig.stage_consistent.sum())}")

if __name__=="__main__": main()
