#!/usr/bin/env python
"""Component plots for Figures 2 to 4 and Supplementary Figures 1 to 3 from the analysis tables."""
from __future__ import annotations
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/"src"))
from biological_staging.config import STAGE_LABELS, T1_DISPLAY, T2_DISPLAY, COL
from biological_staging.io import read_table

DISPLAY_PROTEINS = pd.read_csv(ROOT/"metadata"/"stage_associated_34.csv")["protein"].tolist()


def save(fig,path):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); fig.savefig(path,dpi=300,bbox_inches="tight"); plt.close(fig)


def plot_stage_slopes(table_dir,fig_dir):
    p=Path(table_dir)/"Supplementary_Table_S3A_stage_slopes.csv"
    if not p.exists(): return
    d=pd.read_csv(p).sort_values("stage")
    x=np.linspace(0,3,100); fig,ax=plt.subplots(figsize=(6.5,4.5))
    for _,r in d.iterrows(): ax.plot(x,r.get("intercept_at_time0",0)+r.annual_slope*x,label=r.stage_label)
    ax.set_xlabel("Years from baseline blood draw"); ax.set_ylabel("Model-estimated CDR-SB"); ax.legend(frameon=False,fontsize=8); save(fig,Path(fig_dir)/"Figure3A_CDRSB_trajectory_component.png")


def plot_transition_heatmap(table_dir,fig_dir):
    p=Path(table_dir)/"Figure3B_transition_matrix.csv"
    if not p.exists(): return
    d=pd.read_csv(p); piv=d.pivot(index="baseline_stage",columns="last_stage",values="row_pct")
    fig,ax=plt.subplots(figsize=(6.5,4)); im=ax.imshow(piv.to_numpy(),aspect="auto")
    ax.set_yticks(range(len(piv.index)),[STAGE_LABELS[int(s)] for s in piv.index],fontsize=7); ax.set_xticks(range(len(piv.columns)),[str(c) for c in piv.columns],rotation=45,ha="right",fontsize=7)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]): ax.text(j,i,f"{piv.iloc[i,j]:.0f}%",ha="center",va="center",fontsize=7)
    ax.set_xlabel("Last observed plasma stage"); ax.set_ylabel("Baseline plasma stage"); fig.colorbar(im,ax=ax,label="Row %"); save(fig,Path(fig_dir)/"Figure3B_transition_heatmap_component.png")


def plot_adjacent_heatmap(table_dir,fig_dir):
    p=Path(table_dir)/"Supplementary_Table_S5A_adjacent_136.csv"
    if not p.exists(): return
    d=pd.read_csv(p); d["protein"]=pd.Categorical(d["protein"],categories=DISPLAY_PROTEINS,ordered=True); d=d.sort_values(["protein","stage_from"]); d["contrast"]=d.stage_from.astype(str)+"→"+d.stage_to.astype(str); piv=d.pivot(index="protein",columns="contrast",values="estimate")
    fig,ax=plt.subplots(figsize=(5.5,9)); im=ax.imshow(piv.to_numpy(),aspect="auto")
    ax.set_yticks(range(len(piv.index)),piv.index,fontsize=6); ax.set_xticks(range(len(piv.columns)),piv.columns,fontsize=8); fig.colorbar(im,ax=ax,label="Adjusted difference, SD")
    # mark FDR significant cells
    q=d.pivot(index="protein",columns="contrast",values="q").reindex(index=piv.index,columns=piv.columns)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            if q.iloc[i,j] < .05: ax.text(j,i,"*",ha="center",va="center",fontsize=8)
    save(fig,Path(fig_dir)/"Figure4B_adjacent_stage_heatmap_component.png")


def plot_longitudinal_forest(table_dir,fig_dir):
    p=Path(table_dir)/"Supplementary_Table_S5B_longitudinal_34.csv"
    if not p.exists(): return
    d=pd.read_csv(p); d["protein"]=pd.Categorical(d["protein"],categories=DISPLAY_PROTEINS,ordered=True); d=d.sort_values("protein"); y=np.arange(len(d)); fig,ax=plt.subplots(figsize=(6.5,8))
    ax.errorbar(d.annual_slope_SD_per_year,y,xerr=[d.annual_slope_SD_per_year-d.ci_low,d.ci_high-d.annual_slope_SD_per_year],fmt="o",capsize=2)
    ax.axvline(0,linewidth=.8); ax.set_yticks(y,d.protein,fontsize=6); ax.set_xlabel("Annual protein change, SD/year"); save(fig,Path(fig_dir)/"Figure4D_longitudinal_forest_component.png")


def plot_supp_s3(table_dir,fig_dir):
    p=Path(table_dir)/"Supplementary_Figure_S3_longitudinal_FDR15.csv"
    if not p.exists(): return
    d=pd.read_csv(p); vals=np.column_stack([d.orientation_sign, np.sign(d.annual_slope_SD_per_year)]); fig,ax=plt.subplots(figsize=(4.5,6)); im=ax.imshow(vals,aspect="auto",vmin=-1,vmax=1)
    ax.set_yticks(range(len(d)),d.protein,fontsize=7); ax.set_xticks([0,1],["Stage orientation","Annual change"],rotation=20,ha="right"); save(fig,Path(fig_dir)/"Supplementary_Figure_S3_directional_concordance_component.png")



def plot_t1_auc(table_dir,fig_dir):
    p=Path(table_dir)/"Supplementary_Table_S1A_T1_candidates.csv"
    if not p.exists(): return
    d=pd.read_csv(p).sort_values("AUC")
    y=np.arange(len(d)); fig,ax=plt.subplots(figsize=(6.0,4.3))
    ax.errorbar(d.AUC,y,xerr=[d.AUC-d.AUC_CI_low,d.AUC_CI_high-d.AUC],fmt="o",capsize=2)
    ax.set_yticks(y,d.assay,fontsize=8); ax.set_xlabel("AUC for Aβ PET positivity"); ax.set_xlim(max(0.75,float(d.AUC_CI_low.min())-.02),1.0)
    save(fig,Path(fig_dir)/"Figure2A_T1_AUC_component.png")


def plot_t2_auc(table_dir,fig_dir):
    p=Path(table_dir)/"Supplementary_Table_S2A_T2_candidates.csv"
    if not p.exists(): return
    d=pd.read_csv(p); anchors=list(dict.fromkeys(d.tau_PET_anchor)); markers=list(dict.fromkeys(d.marker))
    fig,ax=plt.subplots(figsize=(6.4,4.2)); x=np.arange(len(anchors))
    offsets=np.linspace(-.18,.18,max(len(markers),1))
    for off,m in zip(offsets,markers):
        g=d[d.marker==m].set_index("tau_PET_anchor").reindex(anchors); ax.plot(x+off,g.AUC,"o-",label=m)
    ax.set_xticks(x,anchors,rotation=15,ha="right",fontsize=8); ax.set_ylabel("AUC"); ax.set_ylim(.65,1.01); ax.legend(frameon=False,fontsize=8)
    save(fig,Path(fig_dir)/"Figure2B_T2_AUC_component.png")


def plot_stage_profile_heatmap(table_dir,fig_dir):
    p=Path(table_dir)/"Figure4A_adjusted_stage_profiles.csv"
    if not p.exists(): return
    d=pd.read_csv(p); cols=[f"adjusted_mean_stage{s}" for s in range(5)]
    # Preserve the manufacturer-domain/protein order from Supplementary Table S4A.
    d["protein"]=pd.Categorical(d["protein"],categories=DISPLAY_PROTEINS,ordered=True); d=d.sort_values("protein"); M=d[cols].to_numpy()
    fig,ax=plt.subplots(figsize=(6.0,8.5)); im=ax.imshow(M,aspect="auto")
    ax.set_yticks(range(len(d)),d.protein,fontsize=6); ax.set_xticks(range(5),[STAGE_LABELS[s] for s in range(5)],rotation=25,ha="right",fontsize=7)
    fig.colorbar(im,ax=ax,label="Age/sex-adjusted standardized mean")
    save(fig,Path(fig_dir)/"Figure4A_stage_profiles_component.png")


def plot_exemplar_profiles(table_dir,fig_dir):
    p=Path(table_dir)/"Figure4A_adjusted_stage_profiles.csv"
    if not p.exists(): return
    d=pd.read_csv(p).set_index("protein"); cols=[f"adjusted_mean_stage{s}" for s in range(5)]; x=np.arange(5)
    for protein in ["GFAP","NfL","pQ-ATXN3"]:
        if protein not in d.index: continue
        fig,ax=plt.subplots(figsize=(4.2,3.2)); ax.plot(x,d.loc[protein,cols].to_numpy(dtype=float),marker="o"); ax.axhline(0,linewidth=.7)
        ax.set_xticks(x,[STAGE_LABELS[s] for s in range(5)],rotation=25,ha="right",fontsize=7); ax.set_ylabel("Adjusted standardized mean"); ax.set_title(protein)
        save(fig,Path(fig_dir)/f"Figure4C_exemplar_{protein.replace('β','beta').replace('/','_')}.png")
    dec=[p for p in ["APOA2","Aβ42","VGF"] if p in d.index]
    if dec:
        centroid=d.loc[dec,cols].mean(axis=0).to_numpy(dtype=float); fig,ax=plt.subplots(figsize=(4.2,3.2))
        for protein in dec: ax.plot(x,d.loc[protein,cols].to_numpy(dtype=float),marker="o",alpha=.55,label=protein)
        ax.plot(x,centroid,marker="o",linewidth=2,label="centroid"); ax.axhline(0,linewidth=.7); ax.set_xticks(x,[STAGE_LABELS[s] for s in range(5)],rotation=25,ha="right",fontsize=7); ax.set_ylabel("Adjusted standardized mean"); ax.legend(frameon=False,fontsize=7)
        save(fig,Path(fig_dir)/"Figure4C_exemplar_decreasing_profiles.png")


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--table-dir",default=str(ROOT/"outputs"/"tables")); ap.add_argument("--figure-dir",default=str(ROOT/"outputs"/"figures")); args=ap.parse_args(); Path(args.figure_dir).mkdir(parents=True,exist_ok=True)
    plot_t1_auc(args.table_dir,args.figure_dir); plot_t2_auc(args.table_dir,args.figure_dir)
    plot_stage_slopes(args.table_dir,args.figure_dir); plot_transition_heatmap(args.table_dir,args.figure_dir)
    plot_stage_profile_heatmap(args.table_dir,args.figure_dir); plot_adjacent_heatmap(args.table_dir,args.figure_dir); plot_exemplar_profiles(args.table_dir,args.figure_dir)
    plot_longitudinal_forest(args.table_dir,args.figure_dir); plot_supp_s3(args.table_dir,args.figure_dir)
    print(f"Figure components written to {args.figure_dir}")
if __name__=="__main__": main()
