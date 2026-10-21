"""Analysis configuration: thresholds, column names and manuscript reference counts."""
from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT_DIR = ROOT / "outputs"
TABLE_DIR = OUTPUT_DIR / "tables"
FIGURE_DIR = OUTPUT_DIR / "figures"
METADATA_DIR = ROOT / "metadata"
REFERENCE_DIR = ROOT / "reference"

# Fixed biomarker thresholds (log2 NPQ; full precision used for computation).
T1_BD_PTAU217 = 9.562758948
T2_EMTBR_EARLY = 9.735217224
T2_EMTBR_MODERATE = 10.35845619
T2_EMTBR_ADVANCED = 11.60847036

# Display values used in manuscript/figures.
T1_DISPLAY = 9.563
T2_DISPLAY = (9.735, 10.358, 11.608)

AMYLOID_CENTILOID_POSITIVE = 24.0
TAU_PET_TEMPORAL_META = 1.376
TAU_PET_NEOC = 1.248
TAU_PET_HIGH_NEOC = 1.6015

STAGE_LABELS = {
    0: "No plasma tau abnormality",
    1: "Tau phosphorylation",
    2: "Early aggregation",
    3: "Moderate aggregation",
    4: "Advanced aggregation",
}
STAGE_ORDER = [0, 1, 2, 3, 4]

# Manuscript reference counts, used only by the QC report.
REFERENCE_COUNTS = {
    "baseline_n": 1035,
    "stageable_n": 1007,
    "unstageable_n": 28,
    "stage_counts": {0: 311, 1: 206, 2: 99, 3: 191, 4: 200},
    "clinical_n": 590,
    "clinical_assessments": 1845,
    "serial_n": 306,
    "serial_draws": 640,
    "baseline_stageable_serial_n": 294,
    "stageable_both_serial_n": 282,
    "transition_progression": 65,
    "transition_stable": 201,
    "transition_regression": 16,
    "progression_eligible_n": 226,
    "progression_eligible_assessments": 627,
    "non_tau_proteins": 197,
    "stage_associated_proteins": 34,
    "longitudinal_fdr_proteins": 15,
    "longitudinal_concordant": 13,
}

# Canonical column names used by the scripts.
COL = {
    "id": "participant_id",
    "draw_id": "draw_id",
    "draw_date": "draw_date",
    "time": "years_from_baseline",
    "age": "age",
    "sex": "sex",
    "education": "education_years",
    "apoe4": "apoe4_carrier",
    "diagnosis": "diagnosis",
    "mmse": "mmse",
    "cdrsb": "cdr_sb",
    "centiloid": "centiloid",
    "amyloid_positive": "amyloid_pet_positive",
    "temporal_meta": "tau_temporal_meta_suvr",
    "neoc": "tau_neoc_suvr",
    "treatment": "anti_amyloid_treatment",
    "phase": "assay_phase",
    "batch": "assay_batch",
    "bd_ptau217": "BD-p-tau217",
    "emtbr_tau": "eMTBR-tau",
}

# Analyte name aliases.
ANALYTE_ALIASES = {
    "NEFL": "NfL",
    "NFL": "NfL",
    "Ab42": "Aβ42",
    "AB42": "Aβ42",
    "pQ_ATXN3": "pQ-ATXN3",
    "pRAB29_T71": "pRAB29-T71",
    "pRAB12_S106": "pRAB12-S106",
    "pPRKN_S65": "pPRKN-S65",
    "mHTT_exon1": "mHTT-exon1",
}
