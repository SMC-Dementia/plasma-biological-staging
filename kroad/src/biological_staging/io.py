"""I/O, schema checks and analyte-name harmonisation."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from .config import ANALYTE_ALIASES


def read_table(path: str | Path, **kwargs) -> pd.DataFrame:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path, **kwargs)
    if suffix in {".tsv"}:
        return pd.read_csv(path, sep="\t", **kwargs)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, **kwargs)
    if suffix == ".parquet":
        return pd.read_parquet(path, **kwargs)
    raise ValueError(f"Unsupported file type: {suffix}")


def require_columns(df: pd.DataFrame, columns: list[str], context: str = "dataset"):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns in {context}: {missing}")


def canonical_analyte(name: str) -> str:
    return ANALYTE_ALIASES.get(str(name), str(name))


def canonicalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    return df.rename(columns={c: canonical_analyte(c) for c in df.columns})


def ensure_output_dir(path: str | Path):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p
