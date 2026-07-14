# Reproduction of Fig. 3 Context-Dependent GRN Construction for Hormone-Dependent Neoplasms

**Bokang Shi**  
Visiting Student, Imperial College London  
Cell Systems Laboratory, Institute for Protein Research, Osaka University  
Summer Internship  
Supervisor: Prof. Mariko Okada  
Date: 30 June 2026

## Introduction

Gene regulatory networks (GRNs) are graphs where genes are represented as nodes and gene-gene relationships are represented as edges. In this paper, the authors compare two related ideas:

- A **context-independent GRN**, which is a general human gene-gene interaction network. It is not specific to any disease.
- A **context-dependent GRN**, which starts from the same general interaction universe but filters and weights interactions according to how relevant the supporting PubMed papers are to a specific disease context.

For this analysis, the disease context was **Neoplasms, Hormone-Dependent**, represented by the MeSH disease identifier `D009376`.

This disease was selected because the broad MeSH parent term `Neoplasms` (`D009369`) was not present as a disease column in the released Fig. 2/Fig. 3 Zenodo data. `Neoplasms, Hormone-Dependent` is an available neoplasm-related MeSH disease term in chunk 1, so it allows a single-disease Fig. 3-style reproduction without changing the underlying GRN construction method.

The aim was to reproduce the Fig. 3-style GRN construction and comparison for hormone-dependent neoplasms: generate a disease-specific context-dependent GRN, compare it against the universal context-independent GRN, and identify central genes using PageRank.

## Summary of what was rerun

The raw Fig. 2 PubMed processing was **not** rerun. In the original workflow, PubMed titles and abstracts are embedded and compared with disease queries using PubMedBERT. That requires the authors' PubMed vector database environment, so this analysis used the released precomputed PubMed-disease similarity scores from Zenodo.

The Fig. 3 GRN construction step **was rerun for Neoplasms, Hormone-Dependent** using an adapted version of the original GitHub script:

Original script:

`context-dependent-GRN/figures/fig3-0_create_context_dependent_GRN.py`

Adapted script:

`context-dependent-GRN/figures/fig3-0_hormone_dependent_neoplasms_create_context_dependent_GRN.py`

The original script is designed to generate context-dependent GRNs for all MeSH diseases. The adapted script keeps the same core Fig. 3 GRN construction logic, but narrows the run to a single disease: Neoplasms, Hormone-Dependent (`D009376`) in chunk 1.

## Data used

The analysis used the authors' released data from Zenodo:

`https://zenodo.org/records/16416117`

The raw BioREX gene interaction file was:

`context-dependent-GRN/data/fig2-3/all_human_gene_interactions_2024-12-18.csv`

This file contains literature-derived human gene-gene interactions. The **context-independent GRN** was generated from all unique gene pairs in this file. This network is disease-agnostic and would be the same regardless of whether the disease being studied was hormone-dependent neoplasms, plasmacytoma, breast cancer, skin cancer, or another disease.

The precomputed PubMedBERT disease-literature similarity file was:

`context-dependent-GRN/data/fig2-3/All_MeSH_diseases_pmid_bert_corrs_chunk1.parquet`

This file contains PubMed article relevance scores for diseases in chunk 1. Neoplasms, Hormone-Dependent is stored in this file as the column `D009376`.

The local MeSH annotation table identifies the disease as:

- disease label: Neoplasms, Hormone-Dependent
- MeSH ID: `D009376`
- chunk number: `1`

## How the GRNs were generated

### Context-independent GRN

The context-independent GRN was generated from:

`all_human_gene_interactions_2024-12-18.csv`

The adapted script took all unique gene-gene pairs and treated each pair as an edge. Each edge was assigned a uniform weight of `1.0`.

This means the context-independent GRN is a universal background network. It is **not disease-specific** and does not change between diseases.

Output file:

`context-dependent-GRN/data/fig2-3/hormone_dependent_neoplasms_fig3_adapted/context_independent_edges.csv`

### Hormone-dependent neoplasms context-dependent GRN

The hormone-dependent neoplasms context-dependent GRN was generated from two inputs:

1. `all_human_gene_interactions_2024-12-18.csv`
2. `All_MeSH_diseases_pmid_bert_corrs_chunk1.parquet`

The adapted script used the `D009376` column from the PubMedBERT score file. For each gene-gene interaction, it looked up the PMID supporting that interaction, mapped the PMID to its hormone-dependent neoplasm relevance score, and retained interactions with relevance score greater than `0.2`.

The retained interactions were then processed using the original Fig. 3 logic:

- remove duplicate PMID-gene-pair entries
- group interactions by gene pair
- sum disease relevance scores per gene pair
- normalize weights to counts per million
- apply `log1p` transformation

Output files:

`context-dependent-GRN/data/fig2-3/hormone_dependent_neoplasms_fig3_adapted/hormone_dependent_neoplasms_context_dependent_edges.csv`

`context-dependent-GRN/data/fig2-3/hormone_dependent_neoplasms_fig3_adapted/hormone_dependent_neoplasms_log1p_CPM_D009376.parquet`

## What the original Fig. 3 script would have done

If run unchanged, `fig3-0_create_context_dependent_GRN.py` would process all MeSH disease chunks:

- chunk 1
- chunk 2
- chunk 3

For every disease column in those chunks, it would generate disease-specific GRN weights. The outputs would be large all-disease matrices:

- `All_MeSH_diseases_log1p_CPM_chunk1.parquet`
- `All_MeSH_diseases_log1p_CPM_chunk2.parquet`
- `All_MeSH_diseases_log1p_CPM_chunk3.parquet`

Each disease would be represented as one column in those matrices. For example:

- `D009376` = Neoplasms, Hormone-Dependent GRN
- `D010954` = plasmacytoma GRN
- other MeSH IDs = other disease-specific GRNs

Therefore, the main purpose of the adaptation was to avoid regenerating every disease and instead run the same Fig. 3 construction logic only for Neoplasms, Hormone-Dependent.

## Code adaptations and justification

The adapted script was:

`context-dependent-GRN/figures/fig3-0_hormone_dependent_neoplasms_create_context_dependent_GRN.py`

The adaptations were:

### 1. Disease restricted to Neoplasms, Hormone-Dependent

What changed: added `DISEASE_ID = "D009376"` and `DISEASE_LABEL = "Neoplasms, Hormone-Dependent"`.

Justification: the broad parent `Neoplasms` (`D009369`) was not present in the released disease score columns, while `D009376` was available as a neoplasm-related MeSH disease term.

### 2. Chunk restricted to chunk 1

What changed: added `CHUNK = 1`.

Justification: the MeSH annotation table and the parquet columns show `D009376` is stored in chunk 1.

The original threshold (`0.2`) and core weighting logic were kept unchanged from `fig3-0`: PMID score mapping, thresholding, duplicate removal, grouping by gene pair, CPM normalization, and `log1p` transformation.

### 3. Context-independent edge list written out

What changed: saved `context_independent_edges.csv`.

Justification: the original script uses the universal pair list internally. Saving it made the baseline GRN explicit for comparison.

### 4. Disease-specific edge list written out

What changed: saved `hormone_dependent_neoplasms_context_dependent_edges.csv`.

Justification: this made the disease-specific GRN directly inspectable.

### 5. PageRank added for both networks

What changed: used `networkx.pagerank(..., weight="weight")`.

Justification: PageRank is used in the Fig. 3 notebooks for representative disease networks; this adaptation applies the same centrality idea to Neoplasms, Hormone-Dependent.

### 6. Overlap table added

What changed: saved `top_pagerank_overlap.csv`.

Justification: this gives a clear comparison between top genes in the two networks.

The original GitHub scripts were not overwritten. The adapted script was added separately so the original repository files remained intact.

## Command run

The adapted Fig. 3 script was run from the extracted GitHub repository's `figures` directory:

```powershell
cd path/to/context-dependent-GRN/figures
python fig3-0_hormone_dependent_neoplasms_create_context_dependent_GRN.py
```

The script generated outputs in:

`context-dependent-GRN/data/fig2-3/hormone_dependent_neoplasms_fig3_adapted/`

## PageRank calculation

PageRank was calculated using NetworkX:

```python
nx.pagerank(graph, weight="weight")
```

For the context-independent GRN, all edge weights were `1.0`, so PageRank reflects centrality in the general background gene-gene interaction network.

For the hormone-dependent neoplasms context-dependent GRN, edge weights were disease-specific literature-derived weights, so PageRank reflects centrality in the hormone-dependent neoplasm-weighted network.

The PageRank outputs were:

- `context_independent_top_pagerank.csv`
- `hormone_dependent_neoplasms_context_dependent_top_pagerank.csv`

## Results

The adapted Fig. 3 script loaded 3,017,208 BioREX interaction rows and identified 891,603 unique gene pairs. After applying the `D009376` relevance threshold, 566,209 interaction rows were retained before aggregation by gene pair.

### Network statistics

| Network | Nodes | Edges | Density | Largest component nodes | Largest component edges |
|---|---:|---:|---:|---:|---:|
| Context-independent GRN | 20,824 | 891,603 | 0.004112 | 20,789 | 891,585 |
| Hormone-dependent neoplasms context-dependent GRN | 15,118 | 242,472 | 0.002122 | 15,051 | 242,438 |

The hormone-dependent neoplasms context-dependent GRN was smaller and less dense than the context-independent GRN. This is expected because disease-context weighting filters the general interaction network to interactions supported by hormone-dependent neoplasm-relevant literature.

### Top PageRank genes

The top 10 PageRank genes in the hormone-dependent neoplasms context-dependent GRN were:

| Rank | Gene | PageRank |
|---:|---|---:|
| 1 | AKT1 | 0.009893 |
| 2 | TP53 | 0.008191 |
| 3 | CTNNB1 | 0.006220 |
| 4 | MAPK1 | 0.006064 |
| 5 | TGFB1 | 0.005761 |
| 6 | MYC | 0.005725 |
| 7 | NFKB1 | 0.005339 |
| 8 | EGFR | 0.004579 |
| 9 | ESR1 | 0.004488 |
| 10 | VEGFA | 0.004188 |

The top 10 PageRank genes in the context-independent GRN were:

| Rank | Gene | PageRank |
|---:|---|---:|
| 1 | AKT1 | 0.004450 |
| 2 | TP53 | 0.003384 |
| 3 | NFKB1 | 0.003179 |
| 4 | MAPK1 | 0.003138 |
| 5 | TNF | 0.003095 |
| 6 | CTNNB1 | 0.002857 |
| 7 | TGFB1 | 0.002856 |
| 8 | MYC | 0.002397 |
| 9 | IL6 | 0.002376 |
| 10 | MTOR | 0.002182 |

### Top gene overlap

Among the top 50 PageRank genes, 42 genes were shared between the context-independent and hormone-dependent neoplasms context-dependent networks.

Genes appearing in the hormone-dependent neoplasms top 50 but not the context-independent top 50 were:

`BRCA1, EREG, EZH2, IGF1, KRAS, PIK3R1, PRL, SNAI1`

Genes appearing in the context-independent top 50 but not the hormone-dependent neoplasms top 50 were:

`APP, CXCL8, FN1, IL10, NFE2L2, PPARG, RELA, TLR4`

## Interpretation

The context-independent GRN is a universal background network. It is identical across disease analyses because it is generated from all unique human gene-gene interactions without disease-specific weighting.

The hormone-dependent neoplasms context-dependent GRN is disease-specific. It starts from the same interaction universe, but interactions are filtered and weighted using disease-specific PubMedBERT relevance scores from the `D009376` column.

Several highly ranked genes were shared between the two networks, including `AKT1`, `TP53`, `MAPK1`, `NFKB1`, `MYC`, `CTNNB1`, and `TGFB1`. These are broad signaling and cancer-associated genes, so it is expected that they remain central in both the general and disease-specific networks.

The disease-specific network also prioritized genes such as `ESR1`, `BRCA1`, `IGF1`, `PRL`, `PIK3R1`, and `SNAI1` within its top 50 genes. This is biologically plausible for hormone-dependent neoplasms because these genes are linked to hormone signaling, cancer growth, or cancer-associated regulatory pathways.

Overall, this analysis reproduced the Fig. 3-style GRN construction and comparison for Neoplasms, Hormone-Dependent by minimally adapting the original GitHub Fig. 3 workflow to a single available neoplasm-related disease context.

## Appendix: exact adapted analysis script

The exact final analysis script added for this project was:

`context-dependent-GRN/figures/fig3-0_hormone_dependent_neoplasms_create_context_dependent_GRN.py`

```python
"""
Hormone-dependent neoplasm-specific adaptation of fig3-0_create_context_dependent_GRN.py.

The original script loops over all MeSH disease columns in all three chunks.
This adaptation keeps the same weighting logic but runs only the MeSH disease
ID D009376, "Neoplasms, Hormone-Dependent", which is stored in chunk 1.
"""

import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


DISEASE_ID = "D009376"
DISEASE_NAME = "hormone_dependent_neoplasms"
DISEASE_LABEL = "Neoplasms, Hormone-Dependent"
CHUNK = 1
THRESHOLD = 0.2
TOP_N = 50


def pair_from_genes(row: pd.Series) -> str:
    return "--".join(sorted([str(row["from_gene"]), str(row["to_gene"])]))


def graph_from_edges(edges: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    graph.add_weighted_edges_from(
        (row.from_gene, row.to_gene, float(row.weight))
        for row in edges.itertuples(index=False)
    )
    return graph


def graph_stats(graph: nx.Graph, label: str) -> dict:
    if graph.number_of_nodes() == 0:
        largest_nodes = 0
        largest_edges = 0
    else:
        largest = graph.subgraph(max(nx.connected_components(graph), key=len))
        largest_nodes = largest.number_of_nodes()
        largest_edges = largest.number_of_edges()

    return {
        "network": label,
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "density": nx.density(graph),
        "largest_component_nodes": largest_nodes,
        "largest_component_edges": largest_edges,
    }


def pagerank_table(graph: nx.Graph, top_n: int) -> pd.DataFrame:
    pagerank = nx.pagerank(graph, weight="weight")
    return (
        pd.DataFrame(pagerank.items(), columns=["gene", "pagerank"])
        .sort_values("pagerank", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )


def main() -> None:
    start = time.time()
    project_root = Path(__file__).resolve().parents[1]
    data_dir = project_root / "data" / "fig2-3"
    output_dir = data_dir / f"{DISEASE_NAME}_fig3_adapted"
    output_dir.mkdir(parents=True, exist_ok=True)

    relations_path = data_dir / "all_human_gene_interactions_2024-12-18.csv"
    pmid_corrs_path = data_dir / f"All_MeSH_diseases_pmid_bert_corrs_chunk{CHUNK}.parquet"

    print(f"Loading BioREX gene interactions: {relations_path}")
    relations_df = pd.read_csv(relations_path)
    relations_df["pair"] = relations_df.apply(pair_from_genes, axis=1)
    all_pairs = np.array(sorted(relations_df["pair"].unique()))
    print("Raw interaction rows:", len(relations_df))
    print("Unique gene pairs:", len(all_pairs))

    print(f"Loading Fig. 2 PubMedBERT scores: {pmid_corrs_path}")
    pmid_corrs = pd.read_parquet(pmid_corrs_path, columns=["pmid", DISEASE_ID])
    pmid_to_corr = pmid_corrs.set_index("pmid")[DISEASE_ID]

    disease_relations_df = relations_df.copy()
    disease_relations_df["corr"] = disease_relations_df["pmid"].map(pmid_to_corr)
    disease_relations_df = disease_relations_df.dropna(subset=["corr"]).reset_index(drop=True)
    disease_relations_df = disease_relations_df[disease_relations_df["corr"] > THRESHOLD].reset_index(drop=True)
    print(f"[{DISEASE_ID}] Interactions after thresholding:", len(disease_relations_df))

    disease_relations_df = disease_relations_df.drop_duplicates(subset=["pmid", "pair"]).reset_index(drop=True)
    disease_edges = (
        disease_relations_df[["pair", "corr"]]
        .groupby("pair")
        .agg({"corr": "sum"})
        .reset_index()
    )
    disease_edges["corr"] = disease_edges["corr"] / disease_edges["corr"].sum() * 1e6
    disease_edges["corr"] = np.log1p(disease_edges["corr"].values)
    disease_edges.columns = ["pair", "weight"]
    disease_edges[["from_gene", "to_gene"]] = disease_edges["pair"].str.split("--", expand=True)
    disease_edges = disease_edges[["from_gene", "to_gene", "weight", "pair"]]

    pair_to_weight = dict(zip(disease_edges["pair"], disease_edges["weight"]))
    matrix_df = pd.DataFrame(
        {
            "pair": all_pairs,
            DISEASE_ID: np.frompyfunc(lambda x: pair_to_weight.get(x, 0), 1, 1)(all_pairs).astype(np.float32),
        }
    )

    context_independent_edges = pd.DataFrame({"pair": all_pairs})
    context_independent_edges[["from_gene", "to_gene"]] = context_independent_edges["pair"].str.split("--", expand=True)
    context_independent_edges["weight"] = 1.0
    context_independent_edges = context_independent_edges[["from_gene", "to_gene", "weight", "pair"]]

    print("Building graphs and calculating PageRank")
    independent_graph = graph_from_edges(context_independent_edges)
    dependent_graph = graph_from_edges(disease_edges)

    stats = pd.DataFrame(
        [
            graph_stats(independent_graph, "context_independent"),
            graph_stats(dependent_graph, f"{DISEASE_NAME}_context_dependent"),
        ]
    )
    independent_pr = pagerank_table(independent_graph, TOP_N)
    dependent_pr = pagerank_table(dependent_graph, TOP_N)

    independent_top = set(independent_pr["gene"])
    dependent_top = set(dependent_pr["gene"])
    overlap = pd.DataFrame(
        {
            "shared_top_genes": [", ".join(sorted(independent_top & dependent_top))],
            f"{DISEASE_NAME}_only_top_genes": [", ".join(sorted(dependent_top - independent_top))],
            "context_independent_only_top_genes": [", ".join(sorted(independent_top - dependent_top))],
            "n_shared": [len(independent_top & dependent_top)],
            "top_n": [TOP_N],
        }
    )

    matrix_df.to_parquet(output_dir / f"{DISEASE_NAME}_log1p_CPM_{DISEASE_ID}.parquet", index=False)
    context_independent_edges.to_csv(output_dir / "context_independent_edges.csv", index=False)
    disease_edges.to_csv(output_dir / f"{DISEASE_NAME}_context_dependent_edges.csv", index=False)
    stats.to_csv(output_dir / "network_stats.csv", index=False)
    independent_pr.to_csv(output_dir / "context_independent_top_pagerank.csv", index=False)
    dependent_pr.to_csv(output_dir / f"{DISEASE_NAME}_context_dependent_top_pagerank.csv", index=False)
    overlap.to_csv(output_dir / "top_pagerank_overlap.csv", index=False)

    elapsed = time.time() - start
    print("\nNetwork stats")
    print(stats.to_string(index=False))
    print(f"\nWrote results to {output_dir}")
    print("Elapsed time:", time.strftime("%H:%M:%S", time.gmtime(elapsed)))


if __name__ == "__main__":
    main()
```
