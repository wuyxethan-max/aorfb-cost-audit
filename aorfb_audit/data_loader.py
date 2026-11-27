"""Data loading and integrity validation for AORFB benchmark datasets."""

from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd

DEFAULT_DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def load_emmel_dataset(data_dir: Path | str | None = None) -> pd.DataFrame:
    """Load the 38-molecule Emmel et al. (2023) audit dataset."""
    root = Path(data_dir) if data_dir is not None else DEFAULT_DATA_DIR
    return pd.read_csv(root / "emmel_38_molecules_audit.csv")


def load_yang_dataset(data_dir: Path | str | None = None) -> pd.DataFrame:
    """Load the Yang et al. (2021) Tables 1–3 hydroquinone cycling dataset."""
    root = Path(data_dir) if data_dir is not None else DEFAULT_DATA_DIR
    return pd.read_csv(root / "yang_tables_1_to_3.csv")


def validate_datasets(
    emmel: pd.DataFrame, yang: pd.DataFrame
) -> Tuple[int, int]:
    """Run structural and physical-bounds assertions on both CSV datasets.

    Returns (n_emmel_molecules, n_yang_rows).
    """
    if yang.groupby("source_table").size().to_dict() != {1: 8, 2: 7, 3: 4}:
        raise ValueError("Yang dataset table counts do not match {1: 8, 2: 7, 3: 4}.")
    if set(yang["molecule"]) != {"O-1", "O-3", "O-5", "asym-O-5"}:
        raise ValueError("Yang dataset molecules do not match expected hydroquinones.")

    if len(emmel) != 38 or not emmel["molecule"].is_unique:
        raise ValueError("Emmel dataset must contain 38 unique molecules.")
    if emmel["electrode_role"].value_counts().to_dict() != {
        "negolyte": 31,
        "posolyte": 7,
    }:
        raise ValueError("Emmel electrode role counts must be 31 negolytes and 7 posolytes.")

    required = [
        "pH_model_fig_S10",
        "potential_fig2_digitized_V",
        "solubility_mol_L",
        "molar_mass_g_mol",
        "n_electrons",
        "fade_pct_day",
        "f_annual_table_S8",
    ]
    if not np.isfinite(emmel[required].to_numpy()).all():
        raise ValueError("Non-finite values found in required Emmel numeric fields.")
    if not (emmel[required[2:]] > 0).all().all():
        raise ValueError("Positive-only molecular properties contain non-positive entries.")
    if not emmel["n_electrons"].isin([1, 2]).all():
        raise ValueError("n_electrons must be 1 or 2 in Emmel dataset.")
    if not emmel["solubility_qualifier"].isin(["=", ">", "~"]).all():
        raise ValueError("Unexpected solubility qualifier in Emmel dataset.")
    if not emmel["fade_qualifier"].isin(["=", "<"]).all():
        raise ValueError("Unexpected fade qualifier in Emmel dataset.")

    eq1_all = (
        emmel["electrolyte_usd_kWh"]
        + emmel["power_usd_kW"] / 4.0
        + emmel["maintenance_usd_kWh"]
    )
    if (eq1_all - emmel["capital_usd_kWh"]).abs().max() >= 0.02:
        raise ValueError("Published Eq. 1 component sum deviates from capital_usd_kWh by >= 0.02.")

    return len(emmel), len(yang)
