#!/usr/bin/env python
"""Generate Table 1 baseline summaries by plasma stage."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL, STAGE_LABELS
from biological_staging.io import read_table, canonicalize_columns
from biological_staging.staging import coerce_stage_value

GROUPS=[("Overall",None)]+[(f"Stage {s} {STAGE_LABELS[s]}",s) for s in range(5)]+[("Unstageable","Unstageable")]

def mean_sd(s):
    x=pd.to_numeric(s,errors="coerce").dropna(); return f"{x.mean():.1f} ({x.std(ddof=1):.1f})" if len(x) else ""
def med_iqr(s):
    x=pd.to_numeric(s,errors="coerce").dropna(); return f"{x.median():.1f} [{x.quantile(.25):.1f}–{x.quantile(.75):.1f}]" if len(x) else ""
def n_pct(mask, denom):
    n=int(mask.sum()); return f"{n} ({100*n/denom:.1f})" if denom else ""

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("output",nargs="?",default=str(ROOT/"outputs"/"tables"/"Table1_baseline.csv")); args=ap.parse_args()
    df=canonicalize_columns(read_table(args.input))
    if "stage" in df.columns:
        df["stage"]=df["stage"].map(coerce_stage_value)
    rows=[]
    specs=[("Age, years",COL["age"],"mean"),("Female sex, n (%)",COL["sex"],"female"),("Education, years",COL["education"],"mean"),("APOE ε4 carrier, n (%)",COL["apoe4"],"binary"),("MMSE score",COL["mmse"],"median"),("CDR-SB",COL["cdrsb"],"median"),("Aβ PET Centiloid",COL["centiloid"],"median"),("Aβ PET positive, n (%)",COL["amyloid_positive"],"binary")]
    for label,col,kind in specs:
        row={"Characteristic":label}
        for gname,g in GROUPS:
            d=df if g is None else df[df.stage==g]
            if col not in d: val=""
            elif kind=="mean": val=mean_sd(d[col])
            elif kind=="median": val=med_iqr(d[col])
            elif kind=="binary": val=n_pct(pd.to_numeric(d[col],errors="coerce")==1, pd.to_numeric(d[col],errors="coerce").notna().sum())
            else:
                female=d[col].astype(str).str.lower().isin(["female","f","1","true"])
                val=n_pct(female,d[col].notna().sum())
            row[gname]=val
        rows.append(row)
    # Diagnosis subrows. Accept common project labels (CU/MCI/Dementia).
    if COL["diagnosis"] in df:
        dx_raw=df[COL["diagnosis"]].astype(str).str.strip().str.lower()
        dx_map={
            "Cognitively unimpaired": {"cu","cognitively unimpaired","cognitively normal","normal"},
            "Mild cognitive impairment": {"mci","mild cognitive impairment"},
            "Dementia": {"dementia","dat","ad dementia","alzheimer dementia","alzheimers dementia"},
        }
        rows.insert(4,{"Characteristic":"Clinical diagnosis, n (%)"})
        for offset,(dx,aliases) in enumerate(dx_map.items(),start=1):
            row={"Characteristic":f"  {dx}"}
            for gname,g in GROUPS:
                d=df if g is None else df[df.stage==g]
                mask=d[COL["diagnosis"]].astype(str).str.strip().str.lower().isin(aliases)
                row[gname]=n_pct(mask,len(d))
            rows.insert(4+offset,row)
    out=pd.DataFrame(rows); Path(args.output).parent.mkdir(parents=True,exist_ok=True); out.to_csv(args.output,index=False); print(f"Wrote {args.output}")
if __name__=="__main__": main()
