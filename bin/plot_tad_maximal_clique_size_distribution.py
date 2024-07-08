#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
# SPDX-License-Identifier: MIT


import argparse
import pathlib
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def make_cli() -> argparse.ArgumentParser:
    def existing_file(arg):
        if (file := pathlib.Path(arg)).exists():
            return file

        raise FileNotFoundError(arg)

    cli = argparse.ArgumentParser(description="Plot maximal clique size distribution across conditions.")

    cli.add_argument(
        "cliques",
        nargs="+",
        type=existing_file,
        help="Path to one or more TSV file with the list of cliques.",
    )
    cli.add_argument(
        "-o",
        "--output-prefix",
        required=True,
        type=pathlib.Path,
        help="Path to output prefix.",
    )
    cli.add_argument(
        "--labels",
        nargs="+",
        type=str,
        help="Sample labels to use for plotting.\n" "When not provided, labels are inferred from input file names.",
    )

    cli.add_argument(
        "-s",
        "--stat",
        type=str,
        choices={"count", "probability", "percent", "density"},
        default="count",
        help="Aggregate statistic to compute in each bin.",
    )

    cli.add_argument(
        "--raise-on-empty-files",
        action="store_true",
        default=False,
        help="Raise an exception when any of the input file(s) is empty.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def save_plot_to_file(fig: plt.Figure, outprefix: pathlib.Path, force: bool, close_after_save: bool = True) -> None:
    png = outprefix.with_suffix(".png")
    svg = outprefix.with_suffix(".svg")
    if not force:
        handle_path_collisions(png, svg)

    outprefix.parent.mkdir(exist_ok=True, parents=True)

    fig.savefig(png, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight")
    if close_after_save:
        plt.close(fig)


def read_cliques(cliques: pathlib.Path, raise_on_empty_files: bool) -> pd.DataFrame:
    df = pd.read_table(cliques).rename(columns={"name": "clique"})
    assert df.columns.tolist() == ["clique", "tad_ids", "size"]

    if raise_on_empty_files and len(df) == 0:
        raise RuntimeError(f"Unable to read any record from {cliques}")

    df["tad_ids"] = df["tad_ids"].str.split(",")

    return df.set_index("clique")


def compute_max_tad_clique_size(cliques: pd.DataFrame) -> pd.DataFrame:
    clique_sizes = {}

    for _, (tads, size) in cliques[["tad_ids", "size"]].iterrows():
        for tad in tads:
            if tad in clique_sizes:
                clique_sizes[tad] = max(clique_sizes[tad], size)
            else:
                clique_sizes[tad] = size

    return pd.DataFrame({"tad": clique_sizes.keys(), "size": clique_sizes.values()})


def plot_maximal_clique_sizes(cliques: Dict[str, pd.DataFrame], stat: str) -> plt.Figure:
    fig, ax = plt.subplots(1, 1)

    data = []
    for label, df in cliques.items():
        data.extend([[label, size] for size in df["size"]])
    df = pd.DataFrame(data, columns=["label", "size"])

    if len(df) > 0:
        sns.histplot(
            df,
            x="size",
            hue="label",
            multiple="dodge",
            ax=ax,
            shrink=0.8,
            discrete=True,
            stat=stat,
            common_norm=False,
        )

        ax.set_xticks(
            range(df["size"].min(), df["size"].max() + 1),
            labels=range(df["size"].min(), df["size"].max() + 1),
        )

    return fig


def generate_labels(paths: Iterable[pathlib.Path]) -> List[str]:
    labels = []
    for p in paths:
        labels.append(str(p.name).rstrip("".join(p.suffixes)))

    return labels


def main():
    args = vars(make_cli().parse_args())

    path_to_cliques = args["cliques"]
    labels = args.get("labels")

    if labels is None:
        labels = generate_labels(path_to_cliques)

    if len(labels) != len(path_to_cliques):
        raise RuntimeError(f"Expected {len(path_to_cliques)} labels, found {len(labels)}")

    cliques = {label: read_cliques(path, args["raise_on_empty_files"]) for label, path in zip(labels, args["cliques"])}

    max_clique_size = {label: compute_max_tad_clique_size(df) for label, df in cliques.items()}

    fig = plot_maximal_clique_sizes(max_clique_size, stat=args["stat"])
    save_plot_to_file(fig, args["output_prefix"], args["force"])


if __name__ == "__main__":
    main()
