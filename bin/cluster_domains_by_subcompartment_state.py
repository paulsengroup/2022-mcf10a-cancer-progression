#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import logging
import pathlib
import sys
from typing import Dict, List, Tuple

import hdbscan
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def make_cli():
    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    def count_or_fraction(arg) -> float:
        n = float(arg)
        if n >= 1 and n.is_integer():
            return n
        if 0 <= n < 1:
            return n

        raise ValueError("Not a count or fraction")

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        nargs="+",
        type=pathlib.Path,
        help="TSV or BED file produced by annotate_domains_with_subcompartments.py.",
    )
    cli.add_argument(
        "--min-cluster-size",
        type=count_or_fraction,
        default=25,
        help="Minimum cluster size.\n"
        "Values in range [0, 1) are interpreted as fractions of the total number of datapoints.\n"
        "See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--min-samples",
        type=count_or_fraction,
        default=5,
        help="Minimum samples.\n"
        "Values in range [0, 1) are interpreted as fractions of the total number of datapoints.\n"
        "See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--cluster-selection-method",
        type=str,
        choices={"eom", "leaf"},
        default="leaf",
        help="See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--cluster-selection-epsilon",
        type=float,
        default=0.1,
        help="See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--distance-metric",
        type=str,
        choices={"euclidean", "manhattan", "chebyshev"},
        default="euclidean",
        help="Distance metric used for clustering.",
    )
    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Path to output prefix (including parent folder(s) but without extension).",
    )
    cli.add_argument(
        "--no-plotting",
        default=False,
        action="store_true",
        help="Skip plot generation. Only output a TSV with the clustered data.",
    )

    cli.add_argument("--plot-title", type=str, default="", help="Title to use for plotting.")

    cli.add_argument("--cmap", type=str, default="deep", help="Colormap used for plotting.")

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


@functools.cache
def get_subcompartment_ranks() -> dict:
    compartment_labels = tuple(["B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3"])
    return {k: v for v, k in enumerate(compartment_labels)}


@functools.cache
def get_compartment_ranks() -> dict:
    compartment_labels = tuple(["B", "A"])
    return {k: v for v, k in enumerate(compartment_labels)}


def import_tsv(path_to_tsvs: List[pathlib.Path]) -> pd.DataFrame:
    stdin = {pathlib.Path("stdin"): sys.stdin, pathlib.Path("-"): sys.stdin}
    return pd.concat((pd.read_table(stdin.get(path, path)) for path in path_to_tsvs))


@functools.cache
def get_color_palette_dict(
    num_colors: int, cmap: str = "deep", outlier_color=(0.5, 0.5, 0.5)
) -> Dict[int, Tuple[float, float, float]]:
    palette = {i: c for i, c in enumerate(sns.color_palette(cmap, n_colors=num_colors))}
    palette[-1] = outlier_color
    return palette


@functools.cache
def get_color_palette_list(
    num_colors: int, cmap: str = "deep", outlier_color=(0.5, 0.5, 0.5)
) -> List[Tuple[float, float, float]]:
    return [
        color
        for _, color in sorted(
            get_color_palette_dict(num_colors, cmap, outlier_color).items(), key=lambda item: item[0]
        )
    ]


def get_cluster_colors(clusters: pd.Series, pvalues: pd.Series, palette: Dict) -> List:
    base_colors = [palette[i] for i in clusters]

    assert len(clusters) == len(pvalues)
    return [sns.desaturate(x, p) for x, p in zip(base_colors, pvalues)]


def plot_scatter(df: pd.DataFrame, suptitle: str, cmap="deep") -> plt.Figure:
    cols = df.filter(regex=r"[AB\d]\.state").columns.tolist()

    color_palette = get_color_palette_dict(num_colors=df["cluster"].max() + 1, cmap=cmap)
    cluster_colors = get_cluster_colors(df["cluster"], df["cluster_pval"], color_palette)

    fig, axs2d = plt.subplots(len(cols), len(cols), figsize=(6.4 * 4, 6.4 * 4))
    for i, (col1, axs1d) in enumerate(zip(cols, axs2d)):
        for j, (col2, ax) in enumerate(zip(cols, axs1d)):
            if j >= i:
                ax.scatter(df[col1], df[col2], alpha=0.2, s=3, c=cluster_colors)
                ax.set(xlim=(-0.1, 1.1), ylim=(-0.1, 1.1))

    legend_patches = [mpl.patches.Patch(facecolor=c, label=f"Cluster #{i}") for i, c in color_palette.items() if i >= 0]
    legend_patches.append(mpl.patches.Patch(facecolor=color_palette[-1], label="Outliers"))

    fig.legend(handles=legend_patches, loc="lower left")
    fig.suptitle(suptitle)
    plt.tight_layout()
    return fig


def rank_clusters(df: pd.DataFrame) -> List[int]:
    """
    Rank clusters based on mean subcompartment rank.
    This can be useful to aid comparison across different plots
    """
    ranked_clusters = {}
    for cluster_id, modal_states in df.groupby("cluster")["state.mode"]:
        ranked_clusters[cluster_id] = modal_states.map(get_subcompartment_ranks()).mean()

    # Rank clusters by value (i.e. by avg state rank)
    ranked_clusters = [int(k) for k, _ in sorted(ranked_clusters.items(), key=lambda item: item[1])]

    # Override order for outlier cluster
    ranked_clusters.remove(-1)
    ranked_clusters.insert(0, -1)

    return ranked_clusters


def plot_distribution(df: pd.DataFrame, suptitle: str, cmap="deep") -> plt.Figure:
    cols = df.filter(regex=r"[AB\d]\.state").columns.tolist()

    num_clusters = df["cluster"].max() + 1
    color_palette = get_color_palette_list(num_colors=num_clusters, cmap=cmap)

    fig, axs = plt.subplots(num_clusters, 1, figsize=(20, 2 * num_clusters))
    for cluster_id, ax, color in zip(rank_clusters(df), axs, color_palette):
        cluster_df = df.loc[df["cluster"] == cluster_id, cols]
        cluster_size = len(cluster_df)
        cluster_size_rel = cluster_size / len(df)
        for _, row in cluster_df.iterrows():
            ax.plot(row, color=color, alpha=0.1)

        if cluster_id >= 0:
            title = f"Cluster #{cluster_id} ({cluster_size} obs; ({100 * cluster_size_rel:.2f}%))"
        else:
            title = f"Outliers ({len(cluster_df)} obs; ({100 * cluster_size_rel:.2f}%))"
        ax.set(title=title, xlim=(0, 1), ylim=(0, 1))

        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels((c.removesuffix(".state") for c in cols))

    fig.suptitle(suptitle)
    plt.tight_layout()
    return fig


def plot_outlier_scores(df: pd.DataFrame, suptitle: str) -> plt.Figure:
    fig, ax = plt.subplots(1, 1)

    scores = df.loc[np.isfinite(df["outlier_score"]), "outlier_score"]
    sns.histplot(scores, kde=True, ax=ax)
    sns.rugplot(scores, ax=ax)

    fig.suptitle(suptitle)
    plt.tight_layout()
    return fig


def run_clustering(
    df: pd.DataFrame, dist_metric, min_cluster_size, min_samples, cluster_selection_method, cluster_selection_epsilon
) -> Tuple[pd.DataFrame, hdbscan.HDBSCAN]:
    cols = df.filter(regex=r"[AB]\d\.state").columns.tolist()
    m = df[cols].to_numpy()
    m = m / m.sum(axis=1)[:, None]

    clusterer = hdbscan.HDBSCAN(
        metric=dist_metric,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method=cluster_selection_method,
        cluster_selection_epsilon=cluster_selection_epsilon,
    )
    clusterer.fit_predict(m)

    df = df.copy()
    df["cluster"] = clusterer.labels_
    df["cluster_pval"] = clusterer.probabilities_
    df["outlier_score"] = clusterer.outlier_scores_
    df[cols] = m

    return df, clusterer


def save_plot_to_file(fig: plt.Figure, outprefix: pathlib.Path, force: bool, close_after_save: bool = True) -> None:
    png = outprefix.with_suffix(".png")
    svg = outprefix.with_suffix(".svg")
    if not force:
        handle_path_collisions(png, svg)

    fig.savefig(png, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight")
    if close_after_save:
        plt.close(fig)


def compute_min_cluster_size(df: pd.DataFrame, n: float) -> int:
    if n >= 1:
        return int(n)

    return int(round(len(df) * n))


def compute_min_samples(df: pd.DataFrame, n: float) -> int:
    return compute_min_cluster_size(df, n)


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def main():
    args = vars(make_cli().parse_args())

    cluster_selection_method = args["cluster_selection_method"]
    cluster_selection_epsilon = args["cluster_selection_epsilon"]

    out_prefix = args["output_prefix"]

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    df = import_tsv(args["tsv"])

    min_cluster_size = compute_min_cluster_size(df, args["min_cluster_size"])
    min_samples = compute_min_samples(df, args["min_samples"])

    df, clusterer = run_clustering(
        df,
        dist_metric=args["distance_metric"],
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method=cluster_selection_method,
        cluster_selection_epsilon=cluster_selection_epsilon,
    )

    outname = pathlib.Path(f"{out_prefix}_clusters.tsv.gz")
    if not args["force"]:
        handle_path_collisions(outname)

    df.to_csv(outname, sep="\t", index=False, header=True)

    print(df["cluster"].unique())

    if not args["no_plotting"]:
        suptitle = "" if args["plot_title"] == "" else args["plot_title"] + " - "
        suptitle += (
            f"HDBSCAN* (min_cluster_size={min_cluster_size}; "
            f"min_samples={min_samples}; "
            f"cluster_selection_method={cluster_selection_method})"
        )

        fig = plot_distribution(df, suptitle=suptitle, cmap=args["cmap"])
        save_plot_to_file(fig, pathlib.Path(f"{out_prefix}_lineplot"), args["force"])

        fig = plot_outlier_scores(df, suptitle=suptitle)
        save_plot_to_file(fig, pathlib.Path(f"{out_prefix}_outliers"), args["force"])

        fig = plot_scatter(df, suptitle=suptitle, cmap=args["cmap"])
        save_plot_to_file(fig, pathlib.Path(f"{out_prefix}_scatter"), args["force"])


if __name__ == "__main__":
    setup_logger()
    main()
