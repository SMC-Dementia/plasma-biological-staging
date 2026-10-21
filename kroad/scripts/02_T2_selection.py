#!/usr/bin/env python
"""T2 marker benchmarking against nested tau-PET anchors (Fig. 2b, Suppl. S2/Table S2)."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from biological_staging.config import COL, TAU_PET_TEMPORAL_META, TAU_PET_NEOC, TAU_PET_HIGH_NEOC
from biological_staging.io import read_table, canonicalize_columns
from biological_staging.stats import auc_ci_bootstrap, bootstrap_youden_cutoff, spearman_with_bootstrap_ci, paired_auc_difference_bootstrap, standardize_within_group
from biological_staging.staging import assign_tau_pet_stage_scalar

MARKERS = ["p-tau205", "BD-p-tau205", "eMTBR-tau"]
ANCHORS = [
    ("Temporal Meta >1.376", COL["temporal_meta"], TAU_PET_TEMPORAL_META),
    ("Neocortical composite >1.248", COL["neoc"], TAU_PET_NEOC),
    ("High neocortical composite >1.6015", COL["neoc"], TAU_PET_HIGH_NEOC),
]


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("input"); ap.add_argument("--outdir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--n-boot-auc",type=int,default=2000,help="Bootstrap replicates for AUC and cutoff (2,000 in the manuscript)"); ap.add_argument("--n-boot-delta",type=int,default=3000,help="Paired delta-AUC bootstrap replicates (3,000 in the manuscript)"); ap.add_argument("--phase-col",default=COL["phase"]); ap.add_argument("--seed",type=int,default=20260826); args=ap.parse_args()
    df=canonicalize_columns(read_table(args.input)); outdir=Path(args.outdir); outdir.mkdir(parents=True,exist_ok=True)
    rows=[]
    for anchor_name, roi, threshold in ANCHORS:
        if roi not in df: continue
        y=(pd.to_numeric(df[roi],errors="coerce")>threshold).astype(float)
        for m in MARKERS:
            if m not in df: continue
            d=pd.DataFrame({"y":y,"x":df[m],"roi":df[roi]}).dropna()
            auc, ci=auc_ci_bootstrap(d.y,d.x,n_boot=args.n_boot_auc,seed=args.seed)
            cut=bootstrap_youden_cutoff(d.y,d.x,n_boot=args.n_boot_auc,seed=args.seed)
            rho,_,rho_ci=spearman_with_bootstrap_ci(d.x,d.roi,n_boot=args.n_boot_auc,seed=args.seed)
            rows.append({"tau_PET_anchor":anchor_name,"positive_N":int(d.y.sum()),"marker":m,"AUC":auc,"AUC_CI_low":ci[0],"AUC_CI_high":ci[1],"Youden_cutoff":cut["cutoff"],"sensitivity":cut["sensitivity"],"specificity":cut["specificity"],"Spearman_rho_vs_ROI":rho,"rho_CI_low":rho_ci[0],"rho_CI_high":rho_ci[1]})
    pd.DataFrame(rows).to_csv(outdir/"Supplementary_Table_S2A_T2_candidates.csv",index=False)

    # Nested PET stage and adjacent-stage AUCs.
    if COL["temporal_meta"] in df and COL["neoc"] in df:
        df["pet_tau_stage"]=[assign_tau_pet_stage_scalar(a,b) for a,b in zip(df[COL["temporal_meta"]],df[COL["neoc"]])]
        order=["Tau PET negative","Early aggregation","Moderate aggregation","Advanced aggregation"]
        out=[]
        for i,s in enumerate(order):
            row={"PET_stage":s,"N":int((df.pet_tau_stage==s).sum())}
            if i>0:
                subset=df[df.pet_tau_stage.isin([order[i-1],s])]
                y=(subset.pet_tau_stage==s).astype(int)
                for m in MARKERS:
                    if m in subset:
                        d=pd.DataFrame({"y":y,"x":subset[m]}).dropna()
                        row[f"{m}_adjacent_AUC"]=roc_auc_score(d.y,d.x)
            out.append(row)
        pd.DataFrame(out).to_csv(outdir/"Supplementary_Table_S2B_PET_hierarchy_adjacent_AUC.csv",index=False)

    # T2 composite sensitivity (Supplementary Fig. 2a): equal-weight mean of
    # within-phase standardised BD-p-tau205 and eMTBR-tau.
    if "eMTBR-tau" in df and "BD-p-tau205" in df:
        z=df.copy()
        cols=["eMTBR-tau","BD-p-tau205"]
        if args.phase_col in z.columns:
            z=standardize_within_group(z,cols,args.phase_col)
        else:
            for c in cols:
                z[c]=(pd.to_numeric(z[c],errors="coerce")-pd.to_numeric(z[c],errors="coerce").mean())/pd.to_numeric(z[c],errors="coerce").std(ddof=1)
        z["T2_equal_weight_composite"]=(z["eMTBR-tau"]+z["BD-p-tau205"])/2
        sens=[]
        for anchor_name,roi,threshold in ANCHORS:
            if roi not in z: continue
            y=(pd.to_numeric(z[roi],errors="coerce")>threshold).astype(float)
            d=pd.DataFrame({"y":y.where(z[roi].notna()),"emtbr":df["eMTBR-tau"],"comp":z["T2_equal_weight_composite"]}).dropna()
            delta,ci,p=paired_auc_difference_bootstrap(d.y,d.comp,d.emtbr,n_boot=args.n_boot_delta,seed=args.seed)
            sens.append({"anchor":anchor_name,"eMTBR_AUC":roc_auc_score(d.y,d.emtbr),"composite_AUC":roc_auc_score(d.y,d.comp),"delta_AUC_composite_minus_eMTBR":delta,"CI_low":ci[0],"CI_high":ci[1],"bootstrap_p":p,"note":"Equal-weight within-phase standardised BD-p-tau205 + eMTBR-tau composite"})
        pd.DataFrame(sens).to_csv(outdir/"Supplementary_Figure_S2A_composite_sensitivity.csv",index=False)
    print(f"T2 outputs written to {outdir}")

if __name__=="__main__": main()
