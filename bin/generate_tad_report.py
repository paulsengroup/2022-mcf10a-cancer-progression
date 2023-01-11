#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib
import warnings
from typing import Dict, List, Tuple

import bioframe as bf
import matplotlib.pyplot as plt
import natsort
import pandas as pd
import seaborn as sns
from matplotlib.patches import Rectangle


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bed",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more BED file(s) with the list of TADs to be processed.",
    )

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Path to output prefix (including parent folder(s) but without extension).",
    )

    cli.add_argument(
        "--labels",
        type=str,
        help="Comma separated list of sample labels to use instead of those inferred from file names.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


def import_bed(path_to_bed: pathlib.Path) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return bf.read_table(path_to_bed, schema="bed")


def import_tads(paths: List[pathlib.Path], labels: List[str]) -> Dict[str, pd.DataFrame]:
    assert len(paths) == len(labels)
    return {label: import_bed(path) for label, path in zip(labels, paths)}


def summarize(df: pd.DataFrame, label: str) -> pd.DataFrame:
    df["size"] = df["end"] - df["start"]

    avg_tad_size = df.groupby("chrom")["size"].mean()
    num_tads = df.groupby("chrom")["size"].count()

    avg_tad_size["all"] = df["size"].mean()
    num_tads["all"] = len(df)

    return pd.DataFrame({f"{label}_avg_tad_size": avg_tad_size, f"{label}_num_tads": num_tads})


def generate_summary(tads: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    df = pd.concat([summarize(df, label) for label, df in tads.items()], axis="columns")
    return (
        df.reset_index().sort_values(by="chrom", key=natsort.natsort_key).fillna(0).round().astype(int, errors="ignore")
    )


def plot_size_distribution(tads: Dict[str, pd.DataFrame], ylim: Tuple[int, int] = (0.0, 3_000_000)) -> plt.Figure:
    labels = []
    sizes = []

    for label, df in tads.items():
        labels.extend([label] * len(df))
        sizes.extend((df["end"] - df["start"]).tolist())

    sizes = pd.DataFrame({"x": labels, "y": sizes})
    sizes["y"] = sizes["y"] / 1.0e6

    counts = {label: len(df) for label, df in tads.items()}

    ylim = [x / 1.0e6 for x in ylim]

    fig, ax = plt.subplots(1, 1)
    sns.violinplot(sizes, x="x", y="y", scale="count", ax=ax, table=counts)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)

    color_palette = sns.color_palette().as_hex()
    artists = [Rectangle((0, 0), 3, 1, facecolor=color_palette[i], edgecolor="black") for i in range(len(counts))]

    ax.legend(
        handles=artists,
        labels=list(counts.values()),
        loc="upper right",
        title="# of TADs",
    )

    ax.set(
        title="TAD size distribution",
        xlabel="",
        ylabel="TAD size (Mbp)",
        ylim=ylim,
    )

    plt.tight_layout()

    return fig


def detect_path_collisions(*args: pathlib.Path) -> None:
    for path in args:
        if path.exists():
            raise RuntimeError(f"Refusing to overwrite file {path}. Pass --force to overwrite existing file(s).")


def main():
    args = vars(make_cli().parse_args())

    output_prefix = args["output_prefix"]
    output_table = output_prefix.with_suffix(".tsv")
    path_to_size_distribution_plot = pathlib.Path(output_prefix.parent) / (
        output_prefix.name + "_size_distribution.png"
    )

    if not args["force"]:
        detect_path_collisions(
            output_table,
            path_to_size_distribution_plot,
            path_to_size_distribution_plot.with_suffix(".svg"),
        )

    labels = args.get("labels")
    paths_to_tads = args["bed"]
    if labels is None:
        labels = [p.name for p in paths_to_tads]
    else:
        labels = labels.split(",")
        if len(labels) != len(paths_to_tads):
            raise RuntimeError(
                f"Mismatch in the number of files and labels: expected {len(paths_to_tads)} labels, found {len(labels)}"
            )

    tads = import_tads(paths_to_tads, labels)

    summary = generate_summary(tads)
    summary.to_csv(output_table, sep="\t", index=False)

    fig = plot_size_distribution(tads)
    fig.savefig(path_to_size_distribution_plot, dpi=600)
    fig.savefig(path_to_size_distribution_plot.with_suffix(".svg"))


if __name__ == "__main__":
    main()
