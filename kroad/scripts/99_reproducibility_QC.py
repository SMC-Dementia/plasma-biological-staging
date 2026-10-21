#!/usr/bin/env python
"""Compare generated outputs with the manuscript reference values in reference/reference_results.json.
Small tolerances are allowed for fitted-model quantities.
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]


def record(rows,item,observed,expected,atol=0.0,rtol=0.0,note=""):
    ok=False
    if observed is not None:
        try:
            ok=bool(np.isclose(float(observed),float(expected),rtol=rtol,atol=atol,equal_nan=False))
        except Exception:
            ok=(observed==expected)
    rows.append({"item":item,"observed":observed,"expected":expected,"pass":ok,"atol":atol,"rtol":rtol,"note":note})


def first_row(df, **conds):
    mask=np.ones(len(df),dtype=bool)
    for col,val in conds.items():
        mask &= (df[col].astype(str)==str(val))
    hit=df.loc[mask]
    return None if hit.empty else hit.iloc[0]


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--table-dir",default=str(ROOT/"outputs"/"tables"))
    ap.add_argument("--output",default=str(ROOT/"outputs"/"QC_REPORT.csv"))
    args=ap.parse_args()
    td=Path(args.table_dir)
    ref=json.loads((ROOT/"reference"/"reference_results.json").read_text(encoding="utf-8"))
    rows=[]
    required=['stage_counts.csv','Supplementary_Table_S1A_T1_candidates.csv','Supplementary_Table_S1B_paired_AUC.csv','Supplementary_Table_S1C_multianalyte.csv','Supplementary_Table_S2A_T2_candidates.csv','Supplementary_Table_S2B_PET_hierarchy_adjacent_AUC.csv','Figure3A_global_stage_by_time.csv','Supplementary_Table_S3A_stage_slopes.csv','Supplementary_Table_S3B_adjacent_slope_contrasts.csv','Figure3B_transition_summary.csv','Figure3C_global_interaction.csv','Figure3C_transition_group_slope_contrasts.csv','Supplementary_Table_S4_full_197.csv','Supplementary_Table_S4_stage_associated_34.csv','Supplementary_Table_S5A_adjacent_136.csv','Supplementary_Table_S5B_longitudinal_34.csv','Supplementary_Table_S5C_domain_scores.csv']
    for filename in required:
        if not (td/filename).exists():
            record(rows,'required_output_'+filename,0,1,note='Required analysis not run or output missing')

    # Stage assignment counts.
    p=td/"stage_counts.csv"
    if p.exists():
        d=pd.read_csv(p)
        for stage,expected in ref["staging"]["stage_counts"].items():
            r=first_row(d,stage=stage)
            record(rows,f"stage_{stage}_count",None if r is None else r["n"],expected)
        r=first_row(d,stage="Unstageable")
        record(rows,"unstageable_count",None if r is None else r["n"],ref["staging"]["unstageable_n"])
        record(rows,"stageable_count",int(d.loc[d.stage.astype(str).isin(["0","1","2","3","4"]),"n"].sum()),ref["staging"]["stageable_n"])
        record(rows,"baseline_count",int(d["n"].sum()),ref["staging"]["baseline_n"])

    # T1 deterministic headline values.
    p=td/"Supplementary_Table_S1A_T1_candidates.csv"
    if p.exists():
        d=pd.read_csv(p)
        for assay,key in [("BD-p-tau217","BD-p-tau217_AUC"),("p-tau217","p-tau217_AUC")]:
            r=first_row(d,assay=assay)
            record(rows,f"T1_{assay}_AUC",None if r is None else r["AUC"],ref["T1"][key],atol=.001)
    p=td/"Supplementary_Table_S1B_paired_AUC.csv"
    if p.exists():
        d=pd.read_csv(p); r=first_row(d,phosphosite="p-tau217")
        record(rows,"T1_p217_delta_AUC",None if r is None else r["delta_AUC_BD_minus_p"],ref["T1"]["delta_AUC_BD_minus_conventional"],atol=.001)
    p=td/"Supplementary_Table_S1C_multianalyte.csv"
    if p.exists():
        d=pd.read_csv(p); r=first_row(d,representation="BD expanded (+212)")
        record(rows,"T1_four_marker_BD_composite_AUC",None if r is None else r["AUC"],ref["T1"]["four_marker_BD_composite_AUC"],atol=.0015)
    p=td/"Supplementary_Table_S1D_T1_threshold.csv"
    if p.exists():
        d=pd.read_csv(p).iloc[0]
        record(rows,"T1_cutoff",d.get("cutoff"),ref["thresholds_full_precision"]["BD-p-tau217_T1"],atol=1e-8)

    # T2 PET anchors and fixed cutoffs.
    p=td/"Supplementary_Table_S2A_T2_candidates.csv"
    if p.exists():
        d=pd.read_csv(p)
        anchor_map={
            "Temporal Meta >1.376":("Temporal Meta",ref["thresholds_full_precision"]["eMTBR_tau_early"]),
            "Neocortical composite >1.248":("NeoC",ref["thresholds_full_precision"]["eMTBR_tau_moderate"]),
            "High neocortical composite >1.6015":("High NeoC",ref["thresholds_full_precision"]["eMTBR_tau_advanced"]),
        }
        for anchor,(rkey,cutoff) in anchor_map.items():
            rr=d[(d["tau_PET_anchor"]==anchor)&(d["marker"]=="eMTBR-tau")]
            if len(rr):
                rr=rr.iloc[0]
                record(rows,f"T2_{rkey}_eMTBR_AUC",rr["AUC"],ref["T2"]["eMTBR_tau_AUCs"][rkey],atol=.001)
                record(rows,f"T2_{rkey}_eMTBR_cutoff",rr["Youden_cutoff"],cutoff,atol=1e-8)
    p=td/"Supplementary_Table_S2B_PET_hierarchy_adjacent_AUC.csv"
    if p.exists():
        d=pd.read_csv(p)
        for label,expected in ref["T2"]["PET_stage_counts"].items():
            rr=d[d["PET_stage"]==label]
            record(rows,f"PET_stage_{label}_n",None if rr.empty else rr.iloc[0]["N"],expected)

    # Longitudinal clinical validation.
    p=td/"Figure3A_global_stage_by_time.csv"
    if p.exists():
        d=pd.read_csv(p).iloc[0]
        record(rows,"clinical_participants",d.get("participants"),ref["clinical"]["participants"])
        record(rows,"clinical_assessments",d.get("assessments"),ref["clinical"]["assessments"])
        record(rows,"clinical_stage_time_chi2",d.get("chi2"),ref["clinical"]["global_stage_time_chi2"],atol=.1)
    p=td/"Supplementary_Table_S3A_stage_slopes.csv"
    if p.exists():
        d=pd.read_csv(p).sort_values("stage")
        for i,e in enumerate(ref["clinical"]["annual_slopes"]):
            rr=d[d.stage==i]
            record(rows,f"stage_{i}_annual_slope",None if rr.empty else rr.iloc[0].annual_slope,e,atol=.005)
    p=td/"Supplementary_Table_S3B_adjacent_slope_contrasts.csv"
    if p.exists():
        d=pd.read_csv(p).sort_values(["stage_from","stage_to"])
        for i,e in enumerate(ref["clinical"]["adjacent_q"]):
            if i<len(d): record(rows,f"adjacent_stage_q_{i}_{i+1}",d.iloc[i].q,e,atol=max(5e-5,abs(e)*.08))

    # Within-person transitions and clinical coupling.
    p=td/"Figure3B_transition_summary.csv"
    if p.exists():
        d=pd.read_csv(p).iloc[0]
        for k,rk in [("baseline_stageable_n","baseline_stageable_serial"),("followup_unstageable_n","followup_unstageable"),("stageable_both_n","stageable_both"),("progression","progression"),("stable","stable"),("regression","regression")]:
            record(rows,k,d[k],ref["transition"][rk])
        record(rows,"directional_forward_fraction",d.get("directional_forward_fraction"),ref["transition"]["forward_fraction_directional"],atol=.001)
        record(rows,"transition_exact_binomial_p",d.get("exact_binomial_p"),ref["transition"]["exact_binomial_p"],atol=1e-9)
    p=td/"Figure3C_global_interaction.csv"
    if p.exists():
        d=pd.read_csv(p).iloc[0]
        record(rows,"progression_eligible_n",d.get("participants"),ref["transition"]["progression_eligible_n"])
        record(rows,"progression_eligible_assessments",d.get("assessments"),ref["transition"]["progression_eligible_assessments"])
    p=td/"Figure3C_transition_group_slope_contrasts.csv"
    if p.exists():
        d=pd.read_csv(p); rr=d[d["contrast"]=="progression vs stable"]
        if len(rr):
            rr=rr.iloc[0]
            record(rows,"progression_vs_stable_slope_difference",rr["estimate"],ref["transition"]["progression_vs_stable_slope_difference"],atol=.006)
            record(rows,"progression_vs_stable_q",rr["q"],ref["transition"]["progression_vs_stable_q"],atol=.004)

    # Cross-sectional / longitudinal proteomics.
    p=td/"Supplementary_Table_S4_stage_associated_34.csv"
    if p.exists(): record(rows,"stage_associated_protein_count",len(pd.read_csv(p)),ref["proteomics"]["stage_associated"])
    p=td/"Supplementary_Table_S5A_adjacent_136.csv"
    if p.exists():
        d=pd.read_csv(p)
        record(rows,"adjacent_contrast_family_size",len(d),ref["proteomics"]["adjacent_contrast_family"])
        first=d[(d.stage_from==0)&(d.stage_to==1)&(d.q<.05)]
        record(rows,"first_boundary_FDR_significant",len(first),ref["proteomics"]["first_boundary_FDR_significant"])
    p=td/"Supplementary_Table_S5B_longitudinal_34.csv"
    if p.exists():
        d=pd.read_csv(p); sig=d[d.q<.05]
        record(rows,"longitudinal_protein_count",len(d),ref["proteomics"]["stage_associated"])
        record(rows,"longitudinal_FDR_count",len(sig),ref["proteomics"]["longitudinal_FDR_significant"])
        record(rows,"longitudinal_concordant_count",int(sig.stage_consistent.sum()),ref["proteomics"]["longitudinal_stage_consistent"])
        if "N_participants" in d:
            record(rows,"longitudinal_participants",int(d.N_participants.max()),ref["proteomics"]["longitudinal_participants"])
        if "N_draws" in d:
            record(rows,"longitudinal_draws",int(d.N_draws.max()),ref["proteomics"]["longitudinal_draws"])
    p=td/"Supplementary_Table_S5C_domain_scores.csv"
    if p.exists():
        d=pd.read_csv(p)
        for domain,sk,qk in [
            ("Inflammation","domain_inflammation_slope","domain_inflammation_q"),
            ("Cytoskeletal abnormalities","domain_cytoskeletal_slope","domain_cytoskeletal_q"),
        ]:
            rr=d[d.domain==domain]
            if len(rr):
                rr=rr.iloc[0]
                record(rows,f"{domain}_slope",rr["annual_slope_SD_per_year"],ref["proteomics"][sk],atol=.004)
                record(rows,f"{domain}_q",rr["q"],ref["proteomics"][qk],atol=max(1e-8,ref["proteomics"][qk]*.25))

    out=pd.DataFrame(rows)
    Path(args.output).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(args.output,index=False)
    if len(out):
        print(out.to_string(index=False))
        print(f"\nQC passed: {int(out['pass'].sum())}/{len(out)} checks")
        if not out["pass"].all():
            print("Some checks did not pass.")
            sys.exit(1)
    else:
        print("No generated outputs found.")
        sys.exit(1)

if __name__=="__main__": main()
