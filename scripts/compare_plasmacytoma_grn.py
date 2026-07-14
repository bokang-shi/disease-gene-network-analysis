import argparse
from pathlib import Path

import networkx as nx
import pandas as pd


DISEASE_ID = "D010954"
DISEASE_NAME = "plasmacytoma"


def find_data_root(repo_root: Path) -> Path:
    candidates = [
        repo_root / "data",
        repo_root / "context-dependent-GRN" / "data",
    ]
    required = Path("fig2-3") / "All_MeSH_diseases_log1p_CPM_chunk1.parquet"
    for candidate in candidates:
        if (candidate / required).exists():
            return candidate
    raise FileNotFoundError(
        "Could not find All_MeSH_diseases_log1p_CPM_chunk1.parquet under data/ "
        "or context-dependent-GRN/data/."
    )


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
    scores = nx.pagerank(graph, weight="weight")
    return (
        pd.DataFrame(scores.items(), columns=["gene", "pagerank"])
        .sort_values("pagerank", ascending=False)
        .head(top_n)
        .reset_index(drop=True)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--top-n", type=int, default=50)
    parser.add_argument("--threshold", type=float, default=0.0)
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    data_root = find_data_root(repo_root)
    out_dir = repo_root / "outputs" / DISEASE_NAME
    out_dir.mkdir(parents=True, exist_ok=True)

    matrix_path = data_root / "fig2-3" / "All_MeSH_diseases_log1p_CPM_chunk1.parquet"
    print(f"Loading {matrix_path}")
    df = pd.read_parquet(matrix_path, columns=["pair", DISEASE_ID])
    df[["from_gene", "to_gene"]] = df["pair"].str.split("--", expand=True)

    context_independent_edges = df[["from_gene", "to_gene"]].copy()
    context_independent_edges["weight"] = 1.0

    context_dependent_edges = (
        df.loc[df[DISEASE_ID] > args.threshold, ["from_gene", "to_gene", DISEASE_ID]]
        .rename(columns={DISEASE_ID: "weight"})
        .copy()
    )

    context_independent_edges.to_csv(
        out_dir / "context_independent_edges.csv", index=False
    )
    context_dependent_edges.to_csv(
        out_dir / "plasmacytoma_context_dependent_edges.csv", index=False
    )

    print("Building graphs")
    independent_graph = graph_from_edges(context_independent_edges)
    dependent_graph = graph_from_edges(context_dependent_edges)

    stats = pd.DataFrame(
        [
            graph_stats(independent_graph, "context_independent"),
            graph_stats(dependent_graph, "plasmacytoma_context_dependent"),
        ]
    )
    stats.to_csv(out_dir / "network_stats.csv", index=False)

    print("Running PageRank")
    independent_pr = pagerank_table(independent_graph, args.top_n)
    dependent_pr = pagerank_table(dependent_graph, args.top_n)
    independent_pr.to_csv(out_dir / "context_independent_top_pagerank.csv", index=False)
    dependent_pr.to_csv(
        out_dir / "plasmacytoma_context_dependent_top_pagerank.csv", index=False
    )

    independent_top = set(independent_pr["gene"])
    dependent_top = set(dependent_pr["gene"])
    overlap = pd.DataFrame(
        {
            "shared_top_genes": [", ".join(sorted(independent_top & dependent_top))],
            "plasmacytoma_only_top_genes": [
                ", ".join(sorted(dependent_top - independent_top))
            ],
            "context_independent_only_top_genes": [
                ", ".join(sorted(independent_top - dependent_top))
            ],
            "n_shared": [len(independent_top & dependent_top)],
            "top_n": [args.top_n],
        }
    )
    overlap.to_csv(out_dir / "top_pagerank_overlap.csv", index=False)

    print("\nNetwork stats")
    print(stats.to_string(index=False))
    print(f"\nWrote results to {out_dir}")


if __name__ == "__main__":
    main()
