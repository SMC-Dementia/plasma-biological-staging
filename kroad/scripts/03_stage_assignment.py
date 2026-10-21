#!/usr/bin/env python
"""Apply the Stage 0–4 algorithm (Fig. 2c; staging backbone)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL, STAGE_LABELS
from biological_staging.io import read_table, canonicalize_columns
from biological_staging.staging import assign_stage

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("output"); ap.add_argument("--bd-col",default=COL["bd_ptau217"]); ap.add_argument("--emtbr-col",default=COL["emtbr_tau"]); args=ap.parse_args()
    df=canonicalize_columns(read_table(args.input)); df=assign_stage(df,args.bd_col,args.emtbr_col)
    df["stage_label"]=df.stage.map(STAGE_LABELS).fillna(df.stage.astype(str))
    out=Path(args.output); out.parent.mkdir(parents=True,exist_ok=True); df.to_csv(out,index=False)
    counts=df["stage"].value_counts(dropna=False)
    rows=[]
    for stage in [0,1,2,3,4,"Unstageable"]:
        rows.append({"stage":stage,"stage_label":STAGE_LABELS.get(stage,"Unstageable"),"n":int(counts.get(stage,0))})
    counts_df=pd.DataFrame(rows)
    counts_df.to_csv(out.parent/"stage_counts.csv",index=False)
    print(counts_df.to_string(index=False))

if __name__=="__main__": main()
