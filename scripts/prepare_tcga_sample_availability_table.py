#!/usr/bin/env python
"""Prepare a compact TCGA tumour/normal sample availability CSV."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUMMARY_CSV = ROOT / "results" / "tcga_degs" / "tcga_deg_summary_compact.csv"
MAPPING_XLSX = ROOT / "data" / "TCGA_MeSH_Mapping_Input_Bokang_Shi.xlsx"
OUTPUT_CSV = ROOT / "results" / "tcga_degs" / "tcga_sample_availability_before_gtex.csv"


def deg_status(n_normal: int) -> str:
    if n_normal >= 2:
        return "Computed"
    return "Skipped"


def main() -> None:
    summary = pd.read_csv(SUMMARY_CSV)
    mapping = pd.read_excel(MAPPING_XLSX, sheet_name="TCGA-MeSH Mapping")

    cancer_names = (
        mapping[["TCGA Abbreviation", "TCGA Disease Name"]]
        .dropna()
        .drop_duplicates(subset=["TCGA Abbreviation"])
        .rename(columns={"TCGA Abbreviation": "TCGA", "TCGA Disease Name": "Cancer Type"})
    )

    table = summary.rename(
        columns={
            "tcga_abbreviation": "TCGA",
            "n_tumor": "Tumour (TCGA)",
            "n_normal": "Normal (TCGA)",
        }
    )[["TCGA", "Tumour (TCGA)", "Normal (TCGA)"]]

    table = table.merge(cancer_names, on="TCGA", how="left")
    table = table[["TCGA", "Cancer Type", "Tumour (TCGA)", "Normal (TCGA)"]]
    table["Tumour (TCGA)"] = pd.to_numeric(table["Tumour (TCGA)"], errors="coerce").fillna(0).astype(int)
    table["Normal (TCGA)"] = pd.to_numeric(table["Normal (TCGA)"], errors="coerce").fillna(0).astype(int)
    table["DEG computation status"] = table["Normal (TCGA)"].map(deg_status)
    table = table.sort_values("TCGA").reset_index(drop=True)

    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUTPUT_CSV, index=False)
    print(OUTPUT_CSV)


if __name__ == "__main__":
    main()
