"""
Extract one context-dependent GRN per selected MeSH ID.

The TCGA-to-MeSH workbook is treated as the extraction manifest. For each unique
MeSH ID marked as available in Zenodo, this script extracts the corresponding
column from the released Fig. 3 all-disease GRN matrix and writes a per-MeSH
edge list. Missing/unavailable MeSH IDs are written to a request table.
"""

from __future__ import annotations

import argparse
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree as ET

import pandas as pd


REQUIRED_COLUMNS = [
    "TCGA Abbreviation",
    "TCGA Disease Name",
    "Mapping Type",
    "Selected MeSH Descriptor",
    "Selected MeSH ID",
    "Available in Zenodo",
    "Fig3 Chunk",
]


def slugify(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_").lower()
    return slug or "unnamed"


def normalize_mesh_id(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"yes", "y", "true", "1"}


def parse_chunk(value: object) -> int | None:
    if pd.isna(value) or str(value).strip() == "":
        return None
    try:
        return int(float(str(value).strip()))
    except ValueError:
        return None


def read_workbook(path: Path, sheet_name: str) -> pd.DataFrame:
    try:
        return pd.read_excel(path, sheet_name=sheet_name)
    except ImportError as exc:
        if "openpyxl" not in str(exc):
            raise
        return read_xlsx_without_openpyxl(path, sheet_name)


def read_xlsx_without_openpyxl(path: Path, sheet_name: str) -> pd.DataFrame:
    """Read a simple XLSX sheet using only the standard library.

    This fallback keeps the script runnable in environments that have pyarrow
    for Parquet files but do not have openpyxl installed.
    """
    ns = {
        "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
        "rel": "http://schemas.openxmlformats.org/package/2006/relationships",
        "office_rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
    }
    with zipfile.ZipFile(path) as zf:
        shared_strings = []
        if "xl/sharedStrings.xml" in zf.namelist():
            root = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in root.findall("main:si", ns):
                text = "".join(t.text or "" for t in si.findall(".//main:t", ns))
                shared_strings.append(text)

        workbook = ET.fromstring(zf.read("xl/workbook.xml"))
        rels = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_targets = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels.findall("rel:Relationship", ns)
        }

        sheet_path = None
        for sheet in workbook.findall("main:sheets/main:sheet", ns):
            if sheet.attrib.get("name") == sheet_name:
                rel_id = sheet.attrib[f"{{{ns['office_rel']}}}id"]
                target = rel_targets[rel_id].lstrip("/")
                sheet_path = Path("xl") / target if not target.startswith("xl/") else Path(target)
                break
        if sheet_path is None:
            available = [s.attrib.get("name", "") for s in workbook.findall("main:sheets/main:sheet", ns)]
            raise ValueError(f"Sheet {sheet_name!r} not found. Available sheets: {available}")

        root = ET.fromstring(zf.read(str(sheet_path).replace("\\", "/")))

    rows: dict[int, dict[int, object]] = defaultdict(dict)
    for cell in root.findall(".//main:sheetData/main:row/main:c", ns):
        ref = cell.attrib.get("r", "")
        match = re.match(r"([A-Z]+)([0-9]+)", ref)
        if not match:
            continue
        col_letters, row_text = match.groups()
        row_idx = int(row_text)
        col_idx = 0
        for ch in col_letters:
            col_idx = col_idx * 26 + ord(ch) - ord("A") + 1

        cell_type = cell.attrib.get("t")
        value = ""
        if cell_type == "inlineStr":
            value = "".join(t.text or "" for t in cell.findall(".//main:t", ns))
        else:
            v = cell.find("main:v", ns)
            if v is not None and v.text is not None:
                if cell_type == "s":
                    value = shared_strings[int(v.text)]
                else:
                    value = v.text
        rows[row_idx][col_idx] = value

    if not rows:
        return pd.DataFrame()

    max_col = max(max(cols) for cols in rows.values())
    table = []
    for row_idx in range(1, max(rows) + 1):
        row = rows.get(row_idx, {})
        table.append([row.get(col_idx, "") for col_idx in range(1, max_col + 1)])

    header = [str(x).strip() for x in table[0]]
    return pd.DataFrame(table[1:], columns=header).replace("", pd.NA)


def combine_unique(values: Iterable[object]) -> str:
    seen = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text and text not in seen:
            seen.append(text)
    return "; ".join(seen)


def build_unique_mesh_table(mapping_df: pd.DataFrame) -> pd.DataFrame:
    missing_columns = [col for col in REQUIRED_COLUMNS if col not in mapping_df.columns]
    if missing_columns:
        raise ValueError(f"Mapping workbook is missing required columns: {missing_columns}")

    df = mapping_df.copy()
    df["Selected MeSH ID"] = df["Selected MeSH ID"].map(normalize_mesh_id)
    df = df[df["Selected MeSH ID"] != ""].copy()
    df["_available"] = df["Available in Zenodo"].map(yes)
    df["_chunk"] = df["Fig3 Chunk"].map(parse_chunk)

    rows = []
    for mesh_id, group in df.groupby("Selected MeSH ID", sort=True):
        descriptor = combine_unique(group["Selected MeSH Descriptor"])
        available = bool(group["_available"].any())
        chunks = sorted({int(x) for x in group["_chunk"].dropna().tolist()})
        rows.append(
            {
                "mesh_id": mesh_id,
                "mesh_descriptor": descriptor,
                "available_in_zenodo": "Yes" if available else "No",
                "fig3_chunk": chunks[0] if chunks else "",
                "tcga_abbreviations": combine_unique(group["TCGA Abbreviation"]),
                "tcga_disease_names": combine_unique(group["TCGA Disease Name"]),
                "mapping_types": combine_unique(group["Mapping Type"]),
                "n_mapping_rows": len(group),
            }
        )
    return pd.DataFrame(rows)


def extract_one_grn(
    mesh_row: pd.Series,
    zenodo_data_dir: Path,
    output_dir: Path,
) -> dict:
    mesh_id = mesh_row["mesh_id"]
    descriptor = mesh_row["mesh_descriptor"]
    chunk = int(mesh_row["fig3_chunk"])
    matrix_path = zenodo_data_dir / f"All_MeSH_diseases_log1p_CPM_chunk{chunk}.parquet"
    if not matrix_path.exists():
        raise FileNotFoundError(f"Missing Zenodo matrix: {matrix_path}")

    df = pd.read_parquet(matrix_path, columns=["pair", mesh_id])
    edges = df.loc[df[mesh_id] > 0, ["pair", mesh_id]].rename(columns={mesh_id: "weight"}).copy()
    edges[["from_gene", "to_gene"]] = edges["pair"].str.split("--", expand=True)
    edges = edges[["from_gene", "to_gene", "weight", "pair"]]

    grn_dir = output_dir / "grns" / f"{mesh_id}__{slugify(descriptor)}"
    grn_dir.mkdir(parents=True, exist_ok=True)
    csv_path = grn_dir / "edges.csv"

    edges.to_csv(csv_path, index=False)

    return {
        "mesh_id": mesh_id,
        "mesh_descriptor": descriptor,
        "status": "extracted",
        "fig3_chunk": chunk,
        "n_edges": int(len(edges)),
        "n_nodes": int(pd.unique(edges[["from_gene", "to_gene"]].values.ravel("K")).size),
        "csv_file": str(csv_path),
        "tcga_abbreviations": mesh_row["tcga_abbreviations"],
        "tcga_disease_names": mesh_row["tcga_disease_names"],
        "mapping_types": mesh_row["mapping_types"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mapping",
        type=Path,
        required=True,
        help="TCGA-to-MeSH mapping workbook.",
    )
    parser.add_argument("--sheet", default="TCGA-MeSH Mapping", help="Workbook sheet name.")
    parser.add_argument(
        "--zenodo-data-dir",
        type=Path,
        default=Path("context-dependent-GRN/data/fig2-3"),
        help="Directory containing All_MeSH_diseases_log1p_CPM_chunk*.parquet.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("outputs/mesh_grns"),
        help="Directory for extracted per-MeSH GRNs and manifests.",
    )
    args = parser.parse_args()

    repo_root = Path.cwd()
    zenodo_data_dir = args.zenodo_data_dir
    if not zenodo_data_dir.is_absolute():
        zenodo_data_dir = repo_root / zenodo_data_dir
    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = repo_root / output_dir

    mapping_df = read_workbook(args.mapping, args.sheet)
    unique_mesh = build_unique_mesh_table(mapping_df)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir = output_dir / "manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    extracted = []
    errors = []
    available_df = unique_mesh[unique_mesh["available_in_zenodo"] == "Yes"].copy()
    for _, row in available_df.iterrows():
        try:
            extracted.append(
                extract_one_grn(
                    row,
                    zenodo_data_dir=zenodo_data_dir,
                    output_dir=output_dir,
                )
            )
            print(f"extracted {row['mesh_id']} {row['mesh_descriptor']}")
        except Exception as exc:  # keep batch moving and report precise failures
            errors.append(
                {
                    "mesh_id": row["mesh_id"],
                    "mesh_descriptor": row["mesh_descriptor"],
                    "status": "error",
                    "error": str(exc),
                }
            )
            print(f"ERROR {row['mesh_id']}: {exc}")

    missing_df = unique_mesh[unique_mesh["available_in_zenodo"] != "Yes"].copy()
    if not missing_df.empty:
        missing_df = missing_df.assign(status="missing_from_zenodo")

    summary_df = pd.DataFrame(extracted + errors)
    summary_df.to_csv(manifest_dir / "extraction_summary.csv", index=False)
    missing_df.to_csv(manifest_dir / "unavailable_mesh_ids.csv", index=False)

    run_metadata = {
        "mapping_workbook": str(args.mapping),
        "sheet": args.sheet,
        "zenodo_data_dir": str(zenodo_data_dir),
        "output_dir": str(output_dir),
        "n_mapping_rows": int(len(mapping_df)),
        "n_unique_mesh_ids": int(len(unique_mesh)),
        "n_available_unique_mesh_ids": int(len(available_df)),
        "n_extracted": int(len(extracted)),
        "n_missing_unique_mesh_ids": int(len(missing_df)),
        "n_errors": int(len(errors)),
    }
    print(run_metadata)


if __name__ == "__main__":
    main()
