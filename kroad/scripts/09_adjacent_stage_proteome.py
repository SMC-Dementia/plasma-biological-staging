#!/usr/bin/env python
"""Age/sex-adjusted adjacent-stage contrasts among the 34 selected proteins (Fig. 4b, Suppl. S5A)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL
from biological_staging.io import read_table, canonicalize_columns
from biological_staging.models import adjacent_stage_contrasts_robust
from biological_staging.stats import bh_fdr
from biological_staging.staging import coerce_stage_value

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("--outdir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--age-col",default=COL["age"]); ap.add_argument("--sex-col",default=COL["sex"]); args=ap.parse_args()
    df=canonicalize_columns(read_table(args.input)); df.stage=df.stage.map(coerce_stage_value); df=df[df.stage.isin(range(5))].copy(); df.stage=df.stage.astype(int)
    meta=read_table(ROOT/"metadata"/"stage_associated_34.csv"); meta["display_order"]=range(len(meta)); proteins=[p for p in meta.protein if p in df.columns]
    rows=[]
    for p in proteins:
        x=pd.to_numeric(df[p],errors="coerce"); z=(x-x.mean())/x.std(ddof=1); d=df[["stage",args.age_col,args.sex_col]].copy(); d[p]=z; d=d.dropna()
        _,cs=adjacent_stage_contrasts_robust(d,p,args.age_col,args.sex_col)
        cs["protein"]=p; rows.append(cs)
    res=pd.concat(rows,ignore_index=True); res["q"]=bh_fdr(res.p); res=res.merge(meta[["protein","domain","display_order"]],on="protein",how="left"); res=res.sort_values(["display_order","stage_from"])
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)
    res.to_csv(outdir/"Supplementary_Table_S5A_adjacent_136.csv",index=False)
    print(f"FDR-significant contrasts: {(res.q<.05).sum()} / {len(res)}")

if __name__=="__main__": main()
