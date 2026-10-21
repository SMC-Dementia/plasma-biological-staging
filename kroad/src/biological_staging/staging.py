"""Plasma stage assignment, tau PET hierarchy and transition labels."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import (
    T1_BD_PTAU217,
    T2_EMTBR_EARLY,
    T2_EMTBR_MODERATE,
    T2_EMTBR_ADVANCED,
    TAU_PET_TEMPORAL_META,
    TAU_PET_NEOC,
    TAU_PET_HIGH_NEOC,
)


def assign_stage_scalar(
    bd_ptau217: float,
    emtbr_tau: float,
    *,
    t1_cutoff: float = T1_BD_PTAU217,
    t2_early: float = T2_EMTBR_EARLY,
    t2_moderate: float = T2_EMTBR_MODERATE,
    t2_advanced: float = T2_EMTBR_ADVANCED,
):
    """Assign Stage 0–4 or ``"Unstageable"``.

    Missing biomarker values return ``np.nan``.
    """
    if pd.isna(bd_ptau217) or pd.isna(emtbr_tau):
        return np.nan
    t1_pos = float(bd_ptau217) >= t1_cutoff
    if not t1_pos:
        return 0 if float(emtbr_tau) < t2_early else "Unstageable"
    if float(emtbr_tau) < t2_early:
        return 1
    if float(emtbr_tau) < t2_moderate:
        return 2
    if float(emtbr_tau) < t2_advanced:
        return 3
    return 4


def assign_stage(df: pd.DataFrame, bd_col: str, emtbr_col: str, out_col: str = "stage") -> pd.DataFrame:
    out = df.copy()
    out[out_col] = [assign_stage_scalar(a, b) for a, b in zip(out[bd_col], out[emtbr_col])]
    return out


def assign_tau_pet_stage_scalar(temporal_meta: float, neoc: float):
    """Nested PET-defined tau hierarchy used for T2 anchoring."""
    if pd.isna(temporal_meta) or pd.isna(neoc):
        return np.nan
    if float(temporal_meta) <= TAU_PET_TEMPORAL_META:
        return "Tau PET negative"
    if float(neoc) <= TAU_PET_NEOC:
        return "Early aggregation"
    if float(neoc) <= TAU_PET_HIGH_NEOC:
        return "Moderate aggregation"
    return "Advanced aggregation"


def transition_label(first_stage, last_stage):
    """Classify first-to-last stage movement."""
    valid = {0, 1, 2, 3, 4}
    if first_stage not in valid:
        return "baseline_unstageable"
    if last_stage not in valid:
        return "followup_unstageable"
    if last_stage > first_stage:
        return "progression"
    if last_stage < first_stage:
        return "regression"
    return "stable"


def coerce_stage_value(x):
    """Coerce CSV-loaded stage values to int 0–4 or 'Unstageable'."""
    if pd.isna(x):
        return np.nan
    s = str(x).strip().lower()
    if s in {"unstageable", "t1-/t2+", "t1−/t2+"}:
        return "Unstageable"
    try:
        value = float(s)
        if not np.isfinite(value) or not value.is_integer():
            return np.nan
        v = int(value)
        return v if v in {0,1,2,3,4} else np.nan
    except Exception:
        return np.nan
