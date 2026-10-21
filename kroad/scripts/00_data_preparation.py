#!/usr/bin/env python
"""Harmonise a baseline table: analyte names, stage assignment and amyloid PET positivity.
Analyte columns are phase/batch-adjusted log2 NPQ values.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from biological_staging.config import AMYLOID_CENTILOID_POSITIVE, COL
from biological_staging.io import read_table, canonicalize_columns, require_columns
from biological_staging.staging import assign_stage


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--id-col", default=COL["id"])
    p.add_argument("--bd-ptau217-col", default=COL["bd_ptau217"])
    p.add_argument("--emtbr-col", default=COL["emtbr_tau"])
    p.add_argument("--centiloid-col", default=COL["centiloid"])
    args = p.parse_args()

    df = canonicalize_columns(read_table(args.input))
    require_columns(df, [args.id_col, args.bd_ptau217_col, args.emtbr_col], "input")
    df = assign_stage(df, args.bd_ptau217_col, args.emtbr_col)
    if args.centiloid_col in df.columns:
        values = pd.to_numeric(df[args.centiloid_col], errors="coerce")
        df[COL["amyloid_positive"]] = (values > AMYLOID_CENTILOID_POSITIVE).where(values.notna()).astype("Int64")
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"Wrote {len(df):,} rows -> {args.output}")

if __name__ == "__main__":
    main()
