"""
Run Fig. 3-style PCA and UMAP on selected MeSH context-dependent GRNs.

This is a minimal script version of figures/fig3-2_all_MeSH_diseases_context_dependent_GRN.ipynb
for the TCGA-selected MeSH GRN matrix already generated in outputs/.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from adjustText import adjust_text
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from umap import UMAP


custom_params = {"axes.spines.right": False, "axes.spines.top": False}
sns.set_theme(style="ticks", rc=custom_params)


def plot_scatter(
    coords: pd.DataFrame,
    x_col: str,
    y_col: str,
    title: str,
    output_path: Path,
    label_points: bool = False,
) -> None:
    fig_size = (8, 7) if label_points else (6, 6)
    fig, ax = plt.subplots(figsize=fig_size)
    sns.scatterplot(
        data=coords,
        x=x_col,
        y=y_col,
        ax=ax,
        color="black",
        s=30,
        alpha=0.8,
    )
    if label_points:
        texts = []
        for row in coords.itertuples(index=False):
            texts.append(
                ax.text(
                    getattr(row, x_col),
                    getattr(row, y_col),
                    row.mesh_id,
                    fontsize=7,
                    alpha=0.9,
                )
            )
        adjust_text(
            texts,
            ax=ax,
            expand=(1.08, 1.25),
            force_text=(0.25, 0.5),
            force_static=(0.2, 0.4),
            force_pull=(0.01, 0.02),
            iter_lim=500,
            arrowprops={"arrowstyle": "-", "color": "0.55", "lw": 0.4},
        )
    ax.set_title(title, loc="left", fontsize=12)
    ax.set_xlabel(x_col)
    ax.set_ylabel(y_col)
    ax.grid(True, linewidth=0.4, alpha=0.35)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--max-pca-components", type=int, default=200)
    args = parser.parse_args()

    matrix = pd.read_parquet(args.matrix)
    if "pair" not in matrix.columns:
        raise ValueError("Input matrix must contain a 'pair' column")

    mesh_ids = [col for col in matrix.columns if col != "pair"]
    relation_corr_df = matrix[mesh_ids].copy()

    is_drop = np.array((relation_corr_df == 0).all(axis=1))
    print("Number of relations to drop:", int(is_drop.sum()))
    relation_corr_df = relation_corr_df[~is_drop].reset_index(drop=True)
    print("Retained relation rows:", len(relation_corr_df))
    print("Disease columns:", len(mesh_ids))
    sparsity = (relation_corr_df.values == 0).mean()
    print(f"Averaged sparsity of the context-dependent GRN: {sparsity:.2%}")

    scaled_corr = StandardScaler().fit_transform(relation_corr_df.values.T)
    n_components = min(args.max_pca_components, scaled_corr.shape[0], scaled_corr.shape[1])
    pca = PCA(n_components=n_components, random_state=args.random_state)
    pca_feat = pca.fit_transform(scaled_corr)
    umap_mapper = UMAP(n_components=2, random_state=args.random_state)
    umap_coord = umap_mapper.fit_transform(pca_feat)

    coords = pd.DataFrame(
        {
            "mesh_id": mesh_ids,
            "PC1": pca_feat[:, 0],
            "PC2": pca_feat[:, 1],
            "UMAP1": umap_coord[:, 0],
            "UMAP2": umap_coord[:, 1],
        }
    )
    explained = pd.DataFrame(
        {
            "component": np.arange(1, n_components + 1),
            "explained_variance_ratio": pca.explained_variance_ratio_,
        }
    )
    run_summary = pd.DataFrame(
        [
            {
                "input_matrix": str(args.matrix),
                "n_mesh_diseases": len(mesh_ids),
                "n_original_relation_rows": len(matrix),
                "n_retained_relation_rows": len(relation_corr_df),
                "n_dropped_all_zero_relation_rows": int(is_drop.sum()),
                "sparsity_after_drop": sparsity,
                "pca_components": n_components,
                "random_state": args.random_state,
            }
        ]
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    coords_path = args.output_dir / "tcga_selected_mesh_grn_pca_umap_coordinates.csv"
    explained_path = args.output_dir / "tcga_selected_mesh_grn_pca_explained_variance.csv"
    summary_path = args.output_dir / "tcga_selected_mesh_grn_pca_umap_summary.csv"
    coords.to_csv(coords_path, index=False)
    explained.to_csv(explained_path, index=False)
    run_summary.to_csv(summary_path, index=False)

    plot_scatter(
        coords,
        "PC1",
        "PC2",
        "PCA of selected MeSH context-dependent GRNs",
        args.output_dir / "tcga_selected_mesh_grn_pca.png",
    )
    plot_scatter(
        coords,
        "PC1",
        "PC2",
        "PCA of selected MeSH context-dependent GRNs",
        args.output_dir / "tcga_selected_mesh_grn_pca_labeled.png",
        label_points=True,
    )
    plot_scatter(
        coords,
        "UMAP1",
        "UMAP2",
        "UMAP of selected MeSH context-dependent GRNs",
        args.output_dir / "tcga_selected_mesh_grn_umap.png",
    )
    plot_scatter(
        coords,
        "UMAP1",
        "UMAP2",
        "UMAP of selected MeSH context-dependent GRNs",
        args.output_dir / "tcga_selected_mesh_grn_umap_labeled.png",
        label_points=True,
    )

    print(f"Wrote coordinates: {coords_path}")
    print(f"Wrote explained variance: {explained_path}")
    print(f"Wrote summary: {summary_path}")


if __name__ == "__main__":
    main()
