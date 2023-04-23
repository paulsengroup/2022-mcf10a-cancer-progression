#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import logging
import pathlib
import pickle
from typing import Dict, List, Tuple, Union

import hdbscan
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.offsetbox import AnchoredText


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        type=pathlib.Path,
        help="Path to *.clusters.tsv.gz produced by cluster_domains_by_subcompartment_state.py.",
    )
    cli.add_argument(
        "pickle",
        type=pathlib.Path,
        help="Path to *.clusterer.pickle.xz produced by cluster_domains_by_subcompartment_state.py.",
    )
    cli.add_argument(
        "--plot-type",
        nargs="+",
        type=str,
        choices={"scatter", "line", "outliers"},
        default="line",
        help="Type(s) of plot to generate.",
    )
    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Path to use as output prefix (including parent folder(s) but without extension).",
    )

    cli.add_argument("--plot-title", type=str, default="HDBSCAN*", help="Title to use for plotting.")
    cli.add_argument("--cmap", type=str, default="deep", help="Colormap used for plotting.")
    cli.add_argument(
        "--label-to-highlight",
        type=str,
        default="all",
        help="Label to higlight. Should be one of the values found in label column.",
    )

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


def import_clusters(path_to_tsv: pathlib.Path, path_to_picke: pathlib.Path) -> Tuple[pd.DataFrame, hdbscan.HDBSCAN]:
    df = pd.read_table(path_to_tsv)
    with open(path_to_picke, "rb") as f:
        return df, pickle.load(f)


@functools.cache
def get_color_palette_dict(
    num_colors: int, cmap: str = "deep", outlier_color=(0.0, 0.0, 0.0)
) -> Dict[int, Tuple[float, float, float]]:
    palette = {i: c for i, c in enumerate(sns.color_palette(cmap, n_colors=num_colors))}
    palette[-1] = outlier_color
    return palette


@functools.cache
def get_color_palette_list(
    num_colors: int, cmap: str = "deep", outlier_color=(0.0, 0.0, 0.0)
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


def plot_distribution(df: pd.DataFrame, suptitle: str, highlight_label: Union[str, None], cmap="deep") -> plt.Figure:
    cols = df.filter(regex=r"[AB\d]\.state").columns.tolist() + ["label"]

    num_clusters = df["cluster"].max() + 1
    color_palette = get_color_palette_list(num_colors=num_clusters, cmap=cmap)

    fig, axs = plt.subplots(num_clusters, 1, figsize=(20, 2 * num_clusters))
    for cluster_id, ax, color in zip(rank_clusters(df), axs, color_palette):
        cluster_df = df.loc[df["cluster"] == cluster_id, cols]
        cluster_size = len(cluster_df)
        cluster_size_rel = cluster_size / len(df)

        cluster_sizes = cluster_df.groupby("label").size()
        ax.add_artist(
            AnchoredText("\n".join((f"{id}: {size}" for id, size in cluster_sizes.items())), loc="upper right")
        )

        if highlight_label is None:
            mask = np.full_like(cluster_df["label"], True, dtype=bool)
        else:
            mask = cluster_df["label"] == highlight_label

        cluster_df.drop(columns="label", inplace=True)

        for _, row in cluster_df[~mask].iterrows():
            ax.plot(row, color="grey", alpha=0.1)

        for _, row in cluster_df[mask].iterrows():
            ax.plot(row, color=color, alpha=0.1)

        if cluster_id >= 0:
            title = f"Cluster #{cluster_id} ({cluster_size} obs; ({100 * cluster_size_rel:.2f}%))"
        else:
            title = f"Outliers ({len(cluster_df)} obs; ({100 * cluster_size_rel:.2f}%))"
        ax.set(title=title, xlim=(0, 1), ylim=(0, 1))

        ax.set_xticks(range(len(cluster_df.columns)))
        ax.set_xticklabels((c.removesuffix(".state") for c in cluster_df.columns))

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

    out_prefix = args["output_prefix"]

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    clusters, clusterer = import_clusters(args["tsv"], args["pickle"])

    dist_metric = clusterer.metric
    min_cluster_size = clusterer.min_cluster_size
    min_samples = clusterer.min_samples
    cluster_selection_method = clusterer.cluster_selection_method
    cluster_selection_epsilon = clusterer.cluster_selection_epsilon

    suptitle = "" if args["plot_title"] == "" else args["plot_title"] + " - "
    suptitle += (
        f"HDBSCAN* (dist_metric={dist_metric}; "
        f"min_cluster_size={min_cluster_size}; "
        f"min_samples={min_samples}; "
        f"cluster_selection_epsilon={cluster_selection_epsilon}; "
        f"cluster_selection_method={cluster_selection_method})"
    )

    plot_types = set(args["plot_type"])
    if "line" in plot_types:
        higlight_label = args["label_to_highlight"]
        if higlight_label.lower() == "all":
            higlight_label = None
        fig = plot_distribution(clusters, suptitle=suptitle, cmap=args["cmap"], highlight_label=higlight_label)
        save_plot_to_file(fig, pathlib.Path(f"{out_prefix}_lineplot"), args["force"])
    if "outliers" in plot_types:
        fig = plot_outlier_scores(clusters, suptitle=suptitle)
        save_plot_to_file(fig, pathlib.Path(f"{out_prefix}_outliers"), args["force"])
    if "scatter" in plot_types:
        fig = plot_scatter(clusters, suptitle=suptitle, cmap=args["cmap"])
        save_plot_to_file(fig, pathlib.Path(f"{out_prefix}_scatter"), args["force"])


if __name__ == "__main__":
    setup_logger()
    main()
