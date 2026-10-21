#!/usr/bin/env python
"""First-to-last within-person stage transition analysis (Fig. 3b)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL, STAGE_LABELS
from biological_staging.io import read_table, require_columns
from biological_staging.staging import coerce_stage_value, transition_label
from biological_staging.stats import exact_directional_binomial

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("--outdir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--id-col",default=COL["id"]); ap.add_argument("--time-col",default=COL["time"]); ap.add_argument("--stage-col",default="stage"); args=ap.parse_args()
    df=read_table(args.input); require_columns(df,[args.id_col,args.time_col,args.stage_col]); df=df.copy(); df[args.stage_col]=df[args.stage_col].map(coerce_stage_value); df=df.sort_values([args.id_col,args.time_col])
    if df.duplicated([args.id_col,args.time_col]).any(): raise ValueError('Duplicate participant and draw time; resolve ordering before transition analysis')
    first=df.drop_duplicates(args.id_col,keep='first')[[args.id_col,args.stage_col]].rename(columns={args.stage_col:"first_stage"})
    last=df.drop_duplicates(args.id_col,keep='last')[[args.id_col,args.stage_col]].rename(columns={args.stage_col:"last_stage"})
    tr=first.merge(last,on=args.id_col); tr=tr[tr.first_stage.isin(range(5))].copy(); tr["transition_group"]=[transition_label(a,b) for a,b in zip(tr.first_stage,tr.last_stage)]
    stageable_both=tr[tr.last_stage.isin(range(5))]
    n_prog=int((stageable_both.transition_group=="progression").sum()); n_stable=int((stageable_both.transition_group=="stable").sum()); n_reg=int((stageable_both.transition_group=="regression").sum()); p=exact_directional_binomial(n_prog,n_reg)

    cols=[0,1,2,3,4,"Unstageable"]
    matrix=[]
    for s in range(5):
        row=tr[tr.first_stage==s]; denom=len(row)
        for dest in cols:
            n=int((row.last_stage==dest).sum()); matrix.append({"baseline_stage":s,"baseline_label":STAGE_LABELS[s],"last_stage":dest,"n":n,"row_pct":100*n/denom if denom else 0})
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)
    tr.to_csv(outdir/"Figure3B_participant_transitions.csv",index=False); pd.DataFrame(matrix).to_csv(outdir/"Figure3B_transition_matrix.csv",index=False)
    pd.DataFrame([{"baseline_stageable_n":len(tr),"followup_unstageable_n":int((tr.transition_group=="followup_unstageable").sum()),"stageable_both_n":len(stageable_both),"progression":n_prog,"stable":n_stable,"regression":n_reg,"directional_forward_fraction":n_prog/(n_prog+n_reg) if n_prog+n_reg else None,"exact_binomial_p":p}]).to_csv(outdir/"Figure3B_transition_summary.csv",index=False)
    print(f"Progression/stable/regression={n_prog}/{n_stable}/{n_reg}; exact binomial P={p:.3g}")

if __name__=="__main__": main()
