import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from biological_staging.staging import assign_stage_scalar, assign_tau_pet_stage_scalar, transition_label
from biological_staging.config import (
    T1_BD_PTAU217,T2_EMTBR_EARLY,T2_EMTBR_MODERATE,T2_EMTBR_ADVANCED,
    TAU_PET_TEMPORAL_META,TAU_PET_NEOC,TAU_PET_HIGH_NEOC,
)


def test_stage_boundaries():
    assert assign_stage_scalar(T1_BD_PTAU217-0.1,T2_EMTBR_EARLY-0.1)==0
    assert assign_stage_scalar(T1_BD_PTAU217-0.1,T2_EMTBR_EARLY)=="Unstageable"
    assert assign_stage_scalar(T1_BD_PTAU217,T2_EMTBR_EARLY-0.1)==1
    assert assign_stage_scalar(T1_BD_PTAU217,T2_EMTBR_EARLY)==2
    assert assign_stage_scalar(T1_BD_PTAU217,T2_EMTBR_MODERATE)==3
    assert assign_stage_scalar(T1_BD_PTAU217,T2_EMTBR_ADVANCED)==4


def test_tau_pet_nested_boundaries():
    assert assign_tau_pet_stage_scalar(TAU_PET_TEMPORAL_META-1e-6, TAU_PET_NEOC-1e-6)=="Tau PET negative"
    assert assign_tau_pet_stage_scalar(TAU_PET_TEMPORAL_META+1e-6, TAU_PET_NEOC-1e-6)=="Early aggregation"
    assert assign_tau_pet_stage_scalar(TAU_PET_TEMPORAL_META+1e-6, TAU_PET_NEOC+1e-6)=="Moderate aggregation"
    assert assign_tau_pet_stage_scalar(TAU_PET_TEMPORAL_META+1e-6, TAU_PET_HIGH_NEOC+1e-6)=="Advanced aggregation"


def test_transition_labels():
    assert transition_label(1,2)=="progression"
    assert transition_label(3,3)=="stable"
    assert transition_label(4,2)=="regression"
    assert transition_label(2,"Unstageable")=="followup_unstageable"


def test_metadata_counts():
    proteins=[x.strip() for x in (ROOT/"metadata"/"non_tau_197.txt").read_text(encoding="utf-8").splitlines() if x.strip()]
    assert len(proteins)==197
    assert len(set(proteins))==197
    import pandas as pd
    selected=pd.read_csv(ROOT/"metadata"/"stage_associated_34.csv")
    assert len(selected)==34
    assert selected["protein"].nunique()==34
