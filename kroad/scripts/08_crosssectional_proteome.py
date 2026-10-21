#!/usr/bin/env python
"""Proteome-wide cross-sectional stage screen (Fig. 4a; Supplementary Table S4)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL
from biological_staging.io import read_table, require_columns, canonicalize_columns
from biological_staging.models import categorical_stage_omnibus, ordinal_stage_robust
from biological_staging.stats import bh_fdr
from biological_staging.staging import coerce_stage_value

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("protein_list",help="One-column CSV or TXT with the 197 eligible non-tau proteins"); ap.add_argument("--outdir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--age-col",default=COL["age"]); ap.add_argument("--sex-col",default=COL["sex"]); args=ap.parse_args()
    df=canonicalize_columns(read_table(args.input)); df["stage"]=df.stage.map(coerce_stage_value); df=df[df.stage.isin(range(5))].copy(); df.stage=df.stage.astype(int)
    plist=read_table(args.protein_list,header=None).iloc[:,0].astype(str).tolist(); proteins=[p for p in plist if p in df.columns]
    if len(proteins)!=len(plist): raise ValueError(f"Incomplete multiplicity family: {len(plist)-len(proteins)} proteins absent from dataset")
    require_columns(df,[args.age_col,args.sex_col,"stage"])
    # Standardise each protein within the baseline stageable cohort.
    for p in proteins: df[p]=(pd.to_numeric(df[p],errors="coerce")-pd.to_numeric(df[p],errors="coerce").mean())/pd.to_numeric(df[p],errors="coerce").std(ddof=1)
    rows=[]
    for p in proteins:
        d=df[[p,"stage",args.age_col,args.sex_col]].dropna()
        fit,om=categorical_stage_omnibus(d,p,args.age_col,args.sex_col)
        ordres=ordinal_stage_robust(d,p,args.age_col,args.sex_col)
        # Age/sex-adjusted stage profile from the categorical model, centred by
        # observed stage frequencies (Figure 4a).
        stage_effects=[0.0] + [float(fit.params.get(f"C(stage)[T.{ss}]", np.nan)) for ss in range(1,5)]
        counts=np.array([(d.stage==ss).sum() for ss in range(5)],dtype=float)
        center=float(np.nansum(counts*np.asarray(stage_effects))/np.nansum(counts))
        adj_means=[v-center for v in stage_effects]
        row={"protein":p,"N":len(d),**om,"ordinal_beta_per_stage_SD":ordres["beta"],"ordinal_SE":ordres["se"],"trend_p":ordres["p"],"Spearman_rho":d[[p,"stage"]].corr(method="spearman").iloc[0,1]}
        row.update({f"adjusted_mean_stage{ss}":adj_means[ss] for ss in range(5)})
        rows.append(row)
    res=pd.DataFrame(rows); res["omnibus_q"]=bh_fdr(res.omnibus_p); res["trend_q"]=bh_fdr(res.trend_p); res=res.sort_values(["omnibus_q","partial_r2"],ascending=[True,False]); res.insert(0,"rank",range(1,len(res)+1))
    domains=read_table(ROOT/"metadata"/"stage_associated_34.csv")[["protein","domain"]].copy(); domains["display_order"]=range(len(domains)); res=res.merge(domains,on="protein",how="left")
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)
    res.to_csv(outdir/"Supplementary_Table_S4_full_197.csv",index=False)
    selected=res[res.omnibus_q<0.05].copy().sort_values(["display_order","rank"],na_position="last")
    selected.to_csv(outdir/"Supplementary_Table_S4_stage_associated_34.csv",index=False)
    profile_cols=["protein","domain","display_order"]+[f"adjusted_mean_stage{ss}" for ss in range(5)]
    selected[profile_cols].to_csv(outdir/"Figure4A_adjusted_stage_profiles.csv",index=False)
    print(f"Tested {len(res)} proteins; q<0.05: {(res.omnibus_q<.05).sum()}")

if __name__=="__main__": main()
