"""
Convert TCGA-selected MeSH PMID-BERT scores into Fig. 3 log1p-CPM GRN matrices.

This is a narrow local adaptation of the authors' figures/fig3-0_create_context_dependent_GRN.py:
the GRN weighting logic is kept the same, while paths and selected MeSH columns
come from the TCGA-MeSH mapping workbook.
"""

from __future__ import annotations

import argparse
import re
import zipfile
from collections import defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

import gc
import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
    "TCGA Abbreviation",
    "TCGA Disease Name",
    "Selected MeSH ID",
    "Available in Zenodo",
    "Fig3 Chunk",
]


def read_workbook(path: Path, sheet_name: str) -> pd.DataFrame:
    try:
        return pd.read_excel(path, sheet_name=sheet_name)
    except ImportError as exc:
        if "openpyxl" not in str(exc):
            raise
        return read_xlsx_without_openpyxl(path, sheet_name)


def read_xlsx_without_openpyxl(path: Path, sheet_name: str) -> pd.DataFrame:
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
                shared_strings.append("".join(t.text or "" for t in si.findall(".//main:t", ns)))

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
                value = shared_strings[int(v.text)] if cell_type == "s" else v.text
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


def normalize_mesh_id(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def availability_is(value: object, expected: str) -> bool:
    return str(value).strip().lower() == expected.lower()


def parse_chunk(value: object) -> int:
    if pd.isna(value) or str(value).strip() == "":
        raise ValueError("Fig3 Chunk is required for Zenodo-available MeSH IDs")
    return int(float(str(value).strip()))


def unique_in_order(values: pd.Series) -> list[str]:
    seen = set()
    out = []
    for value in values:
        mesh_id = normalize_mesh_id(value)
        if mesh_id and mesh_id not in seen:
            seen.add(mesh_id)
            out.append(mesh_id)
    return out


def selected_mapping(mapping: pd.DataFrame, availability: str) -> pd.DataFrame:
    missing = [col for col in REQUIRED_COLUMNS if col not in mapping.columns]
    if missing:
        raise ValueError(f"Mapping workbook is missing columns: {missing}")
    df = mapping[mapping["Available in Zenodo"].map(lambda x: availability_is(x, availability))].copy()
    df["Selected MeSH ID"] = df["Selected MeSH ID"].map(normalize_mesh_id)
    df = df[df["Selected MeSH ID"] != ""].copy()
    return df


def build_all_pairs(relations_df: pd.DataFrame) -> np.ndarray:
    return np.array(
        sorted(
            list(
                set(
                    relations_df.apply(
                        lambda x: "--".join(sorted([x["from_gene"], x["to_gene"]])),
                        axis=1,
                    )
                )
            )
        )
    )


def append_fig3_weights(
    chunk_results: dict[str, object],
    summary_rows: list[dict[str, object]],
    relations_df: pd.DataFrame,
    all_pairs: np.ndarray,
    pmid_corrs: pd.DataFrame,
    disease_columns: list[str],
    source_file: Path,
    threshold: float,
) -> None:
    available = set(pmid_corrs.columns)
    missing = [col for col in disease_columns if col not in available]
    if missing:
        raise ValueError(f"{source_file} is missing requested MeSH columns: {missing}")

    for col in disease_columns:
        disease_relations_df = relations_df.copy()
        disease_relations_df["corr"] = disease_relations_df["pmid"].map(pmid_corrs.set_index("pmid")[col])
        disease_relations_df = disease_relations_df.dropna(subset=["corr"]).reset_index(drop=True)
        disease_relations_df = disease_relations_df[disease_relations_df["corr"] > threshold].reset_index(drop=True)
        print(f"[{col}] Number of interactions after thresholding:", len(disease_relations_df))

        n_relations_after_threshold = len(disease_relations_df)
        disease_relations_df["pair"] = disease_relations_df.apply(
            lambda x: "--".join(sorted([x["from_gene"], x["to_gene"]])),
            axis=1,
        )
        disease_relations_df = disease_relations_df.drop_duplicates(subset=["pmid", "pair"]).reset_index(drop=True)

        disease_relations_df = (
            disease_relations_df[["pair", "corr"]]
            .groupby("pair")
            .agg({"corr": "sum"})
            .reset_index()
        )
        disease_relations_df["corr"] = disease_relations_df["corr"] / disease_relations_df["corr"].sum() * 1e6
        disease_relations_df["corr"] = np.log1p(disease_relations_df["corr"].values)
        disease_relations_df.columns = ["pair", "weight"]
        pair2weight = dict(zip(disease_relations_df["pair"], disease_relations_df["weight"]))
        weights = np.frompyfunc(lambda x: pair2weight.get(x, 0), 1, 1)(all_pairs).astype(np.float32)
        chunk_results[col] = weights

        summary_rows.append(
            {
                "mesh_id": col,
                "source_file": str(source_file),
                "threshold": threshold,
                "n_interactions_after_threshold": n_relations_after_threshold,
                "n_unique_weighted_pairs": int((weights > 0).sum()),
            }
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--sheet", default="TCGA-MeSH Mapping")
    parser.add_argument("--availability", choices=["Yes", "No"], required=True)
    parser.add_argument(
        "--relations",
        type=Path,
        default=Path("context-dependent-GRN/data/fig2-3/all_human_gene_interactions_2024-12-18.csv"),
    )
    parser.add_argument(
        "--all-mesh-pmid-dir",
        type=Path,
        default=Path("context-dependent-GRN/data/fig2-3"),
    )
    parser.add_argument("--missing-pmid-bert", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--threshold", type=float, default=0.2)
    args = parser.parse_args()

    mapping = read_workbook(args.mapping, args.sheet)
    selected = selected_mapping(mapping, args.availability)
    if selected.empty:
        raise ValueError(f"No mapping rows found for Available in Zenodo = {args.availability}")

    print(f"Loading BioREX gene interactions: {args.relations}")
    relations_df = pd.read_csv(args.relations)
    all_pairs = build_all_pairs(relations_df)
    print("Raw interaction rows:", len(relations_df))
    print("Unique gene pairs:", len(all_pairs))

    chunk_results: dict[str, object] = {"pair": all_pairs}
    summary_rows: list[dict[str, object]] = []

    if args.availability == "Yes":
        selected["_chunk"] = selected["Fig3 Chunk"].map(parse_chunk)
        for chunk, group in selected.groupby("_chunk", sort=True):
            disease_columns = unique_in_order(group["Selected MeSH ID"])
            source_file = args.all_mesh_pmid_dir / f"All_MeSH_diseases_pmid_bert_corrs_chunk{chunk}.parquet"
            print(f"Loading Fig. 2 PubMedBERT scores: {source_file}")
            pmid_corrs = pd.read_parquet(source_file, columns=["pmid", *disease_columns])
            append_fig3_weights(
                chunk_results,
                summary_rows,
                relations_df,
                all_pairs,
                pmid_corrs,
                disease_columns,
                source_file,
                args.threshold,
            )
            del pmid_corrs
            gc.collect()
    else:
        if args.missing_pmid_bert is None:
            raise ValueError("--missing-pmid-bert is required when --availability No")
        disease_columns = unique_in_order(selected["Selected MeSH ID"])
        print(f"Loading Fig. 2 PubMedBERT scores: {args.missing_pmid_bert}")
        pmid_corrs = pd.read_parquet(args.missing_pmid_bert, columns=["pmid", *disease_columns])
        append_fig3_weights(
            chunk_results,
            summary_rows,
            relations_df,
            all_pairs,
            pmid_corrs,
            disease_columns,
            args.missing_pmid_bert,
            args.threshold,
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(chunk_results).to_parquet(args.output, index=False)
    pd.DataFrame(summary_rows).to_csv(args.summary_output, index=False)
    print(f"Wrote log1p-CPM matrix: {args.output}")
    print(f"Wrote summary: {args.summary_output}")


if __name__ == "__main__":
    main()
