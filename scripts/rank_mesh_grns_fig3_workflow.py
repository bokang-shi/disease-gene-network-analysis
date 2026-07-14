"""
Rank genes for MeSH context-dependent GRNs using the authors' Fig. 3 workflow.

This follows figures/fig3-1_repr_diseases_context_dependent_GRN.ipynb:
1. Build a weighted graph from one disease-specific GRN column.
2. Extract the largest connected component.
3. Apply degree penalization with beta=1.0:
   W = D^-beta (C + A) D^-beta, where C is common-neighbor counts and A is
   the unweighted adjacency matrix.
4. Multiply W by the weighted adjacency and run NetworkX PageRank.

The notebook materializes dense matrices. This script performs the same matrix
operations sparsely so it can run across all selected MeSH diseases.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp


def split_pairs(pair_series: pd.Series) -> tuple[np.ndarray, np.ndarray]:
    pairs = pair_series.str.split("--", n=1, expand=True)
    if pairs.shape[1] != 2:
        raise ValueError("Expected pair values in 'geneA--geneB' format")
    return pairs[0].to_numpy(), pairs[1].to_numpy()


def graph_from_column(
    from_genes: np.ndarray,
    to_genes: np.ndarray,
    weights: np.ndarray,
) -> nx.Graph:
    keep = weights > 0
    graph = nx.Graph()
    graph.add_weighted_edges_from(
        (u, v, float(w))
        for u, v, w in zip(from_genes[keep], to_genes[keep], weights[keep])
    )
    return graph


def largest_component(graph: nx.Graph) -> nx.Graph:
    connected_components = list(nx.connected_components(graph))
    if not connected_components:
        return graph.copy()
    return graph.subgraph(connected_components[np.argmax([len(c) for c in connected_components])]).copy()


def degree_penalized_graph(graph: nx.Graph, beta: float) -> nx.Graph:
    nodes = list(graph.nodes)
    if not nodes:
        return graph.copy()

    adjacency = nx.to_scipy_sparse_array(graph, nodelist=nodes, weight=None, format="csr", dtype=np.float64)
    weighted_adjacency = nx.to_scipy_sparse_array(
        graph,
        nodelist=nodes,
        weight="weight",
        format="csr",
        dtype=np.float64,
    )

    common_neighbors = adjacency.T @ adjacency
    common_neighbors = common_neighbors - sp.diags(common_neighbors.diagonal(), format="csr")
    c_prime = common_neighbors + adjacency

    degrees = np.asarray(adjacency.sum(axis=1)).ravel()
    d_inv_beta = sp.diags(np.power(degrees, -beta), format="csr")
    penalty = d_inv_beta.T @ c_prime @ d_inv_beta
    degree_penalized_adjacency = penalty.multiply(weighted_adjacency).tocsr()

    penalized_graph = nx.from_scipy_sparse_array(degree_penalized_adjacency, create_using=nx.Graph)
    return nx.relabel_nodes(penalized_graph, dict(enumerate(nodes)))


def pagerank_table(graph: nx.Graph, mesh_id: str) -> pd.DataFrame:
    pagerank = nx.pagerank(graph)
    out = (
        pd.DataFrame(pagerank.items(), columns=["gene", "pagerank"])
        .sort_values("pagerank", ascending=False)
        .reset_index(drop=True)
    )
    out.insert(0, "mesh_id", mesh_id)
    out.insert(1, "rank", np.arange(1, len(out) + 1))
    return out


def rank_one_mesh(
    mesh_id: str,
    from_genes: np.ndarray,
    to_genes: np.ndarray,
    weights: np.ndarray,
    beta: float,
) -> tuple[pd.DataFrame, dict[str, object]]:
    start = time.time()
    graph = graph_from_column(from_genes, to_genes, weights)
    main_subgraph = largest_component(graph)
    penalized_graph = degree_penalized_graph(main_subgraph, beta)
    ranking = pagerank_table(penalized_graph, mesh_id)
    elapsed = time.time() - start

    stats = {
        "mesh_id": mesh_id,
        "raw_nodes": graph.number_of_nodes(),
        "raw_edges": graph.number_of_edges(),
        "raw_density": nx.density(graph) if graph.number_of_nodes() > 1 else 0.0,
        "largest_component_nodes": main_subgraph.number_of_nodes(),
        "largest_component_edges": main_subgraph.number_of_edges(),
        "ranked_genes": len(ranking),
        "beta": beta,
        "elapsed_seconds": elapsed,
    }
    return ranking, stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-rankings", type=Path, required=True)
    parser.add_argument("--output-stats", type=Path, required=True)
    parser.add_argument("--top-n-csv", type=Path)
    parser.add_argument("--beta", type=float, default=1.0)
    args = parser.parse_args()

    matrix = pd.read_parquet(args.matrix)
    if "pair" not in matrix.columns:
        raise ValueError("Input matrix must contain a 'pair' column")

    mesh_ids = [col for col in matrix.columns if col != "pair"]
    from_genes, to_genes = split_pairs(matrix["pair"])

    rankings = []
    stats_rows = []
    for mesh_id in mesh_ids:
        print(f"Ranking {mesh_id}")
        ranking, stats = rank_one_mesh(
            mesh_id,
            from_genes,
            to_genes,
            matrix[mesh_id].to_numpy(),
            args.beta,
        )
        rankings.append(ranking)
        stats_rows.append(stats)
        print(
            f"[{mesh_id}] ranked_genes={stats['ranked_genes']} "
            f"largest_edges={stats['largest_component_edges']} "
            f"elapsed={stats['elapsed_seconds']:.1f}s"
        )

    all_rankings = pd.concat(rankings, ignore_index=True)
    stats_df = pd.DataFrame(stats_rows)

    args.output_rankings.parent.mkdir(parents=True, exist_ok=True)
    args.output_stats.parent.mkdir(parents=True, exist_ok=True)
    all_rankings.to_parquet(args.output_rankings, index=False)
    stats_df.to_csv(args.output_stats, index=False)
    if args.top_n_csv is not None:
        args.top_n_csv.parent.mkdir(parents=True, exist_ok=True)
        top_n = all_rankings[all_rankings["rank"] <= 100].copy()
        top_n.to_csv(args.top_n_csv, index=False)

    print(f"Wrote rankings: {args.output_rankings}")
    print(f"Wrote stats: {args.output_stats}")
    if args.top_n_csv is not None:
        print(f"Wrote top-100 CSV: {args.top_n_csv}")


if __name__ == "__main__":
    main()
