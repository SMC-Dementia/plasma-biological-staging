#!/usr/bin/env python
"""Clinical coupling of stage movement (Fig. 3c)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import COL
from biological_staging.io import read_table, require_columns
from biological_staging.models import fit_mixedlm, fixed_effect_params_and_cov, likelihood_ratio_test
from biological_staging.stats import linear_contrast, bh_fdr
from biological_staging.staging import coerce_stage_value

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("clinical_long"); ap.add_argument("transitions"); ap.add_argument("--outdir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--id-col",default=COL["id"]); ap.add_argument("--time-col",default=COL["time"]); ap.add_argument("--cdrsb-col",default=COL["cdrsb"]); args=ap.parse_args()
    c=read_table(args.clinical_long); t=read_table(args.transitions); require_columns(t,[args.id_col,"first_stage","last_stage","transition_group"])
    t["first_stage"]=t.first_stage.map(coerce_stage_value); t=t[t.first_stage.isin([0,1,2,3]) & t.transition_group.isin(["progression","stable","regression"])].copy()
    d=c.merge(t[[args.id_col,"first_stage","transition_group"]],on=args.id_col,how="inner")
    d=d.dropna(subset=[args.id_col,args.time_col,args.cdrsb_col,"first_stage","transition_group"]).copy()
    d["first_stage"]=d["first_stage"].astype(int)
    d["transition_group"]=pd.Categorical(d.transition_group,categories=["stable","progression","regression"])
    formula=f"Q('{args.cdrsb_col}') ~ {args.time_col} * C(transition_group) + {args.time_col} * C(first_stage)"
    reduced=f"Q('{args.cdrsb_col}') ~ {args.time_col} + C(transition_group) + {args.time_col} * C(first_stage)"
    fit=fit_mixedlm(formula,d,args.id_col,re_formula=f"~{args.time_col}",reml=False); red=fit_mixedlm(reduced,d,args.id_col,re_formula=f"~{args.time_col}",reml=False)
    lr,dfd,p_global=likelihood_ratio_test(fit,red,2); params,cov=fixed_effect_params_and_cov(fit)
    rows=[]
    for grp in ["progression","regression"]:
        term=f"{args.time_col}:C(transition_group)[T.{grp}]"; cst=linear_contrast(params,cov,{term:1.0}); cst.update({"contrast":f"{grp} vs stable"}); rows.append(cst)
    out=pd.DataFrame(rows); out["q"]=bh_fdr(out.p)
    outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True); out.to_csv(outdir/"Figure3C_transition_group_slope_contrasts.csv",index=False); pd.DataFrame([{"chi2":lr,"df":dfd,"p":p_global,"participants":d[args.id_col].nunique(),"assessments":len(d)}]).to_csv(outdir/"Figure3C_global_interaction.csv",index=False)
    print(out[["contrast","estimate","ci_low","ci_high","q"]].to_string(index=False))

if __name__=="__main__": main()
