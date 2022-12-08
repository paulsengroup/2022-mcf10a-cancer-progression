#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import itertools
import pathlib
from typing import List, Union

import bioframe as bf
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bedgraph",
        type=pathlib.Path,
        help="Path to the (sub)compartment bedGraph generated by dcHiC.",
    )

    cli.add_argument(
        "output-prefix",
        type=pathlib.Path,
        help="Path to output prefix (including parent folder(s) but without extension).",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


def import_data(path_to_bedgraph: pathlib.Path) -> pd.DataFrame:
    df = pd.read_table(path_to_bedgraph).rename(columns={"chr": "chrom"})
    df1 = df[["chrom", "start", "end", "padj"]].copy()
    df2 = df.filter(regex=".state$")

    df1[df2.columns] = df2
    return df1


@functools.cache
def get_compartment_ranks() -> dict:
    compartment_labels = tuple(
        ["B", "B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3", "A"]
    )
    return {k: v for v, k in enumerate(compartment_labels)}


def get_compartment_labels_from_df(df: pd.DataFrame) -> List[str]:
    labels = (
        pd.concat([df[col] for col in df.columns if col.endswith(".state")])
        .unique()
        .tolist()
    )

    return list(
        sorted(labels, key=lambda x: [get_compartment_ranks().get(y) for y in x])
    )


def plot_compartment_size_distribution_by_condition(
    df: pd.DataFrame, axs: List[plt.Axes], filter: Union[str, None] = None
) -> None:
    if filter is not None:
        df = df[df["chrom"] == filter]

    cols = [col for col in df.columns if col.endswith(".state")]
    assert len(cols) == len(axs)

    for ax, col in zip(axs, cols):
        df1 = bf.merge(df[["chrom", "start", "end", col]], min_dist=0, on=[col]).drop(
            columns=["n_intervals"]
        )

        df1["Subcompartment length"] = df1["end"] - df1["start"]

        sns.violinplot(df1, x=col, y="Subcompartment length", ax=ax)
        sns.violinplot(df1, x=col, y="Subcompartment length", ax=ax)
        ax.set(title=col, ylim=[0, int(5.0e6)])


def plot_compartment_size_distribution_by_subcomp(
    df: pd.DataFrame, axs: List[plt.Axes], filter: Union[str, None] = None
) -> None:
    subcomps = get_compartment_labels_from_df(df)
    assert len(subcomps) == len(axs)

    if filter is not None:
        df = df[df["chrom"] == filter]

    cols = [col for col in df.columns if col.endswith(".state")]

    dfs = []
    for col in cols:
        df1 = (
            bf.merge(df[["chrom", "start", "end", col]], min_dist=0, on=[col])
            .drop(columns=["n_intervals"])
            .rename(columns={col: "Subcompartment"})
        )

        df1["Condition"] = col
        dfs.append(df1)

    df = pd.concat(dfs)
    df["Subcompartment length"] = df["end"] - df["start"]

    for subcmp, ax in zip(subcomps, axs):
        sns.violinplot(
            df[df["Subcompartment"] == subcmp],
            x="Condition",
            y="Subcompartment length",
            ax=ax,
        )
        ax.set(title=subcmp, ylim=[0, int(5.0e6)])


def main():
    args = vars(make_cli().parse_args())

    df = import_data(args["bedgraph"])
    outprefix = args["output-prefix"]

    fig, axs = plt.subplots(2, 4, figsize=(4 * 6.4, 2 * 6.4))
    plot_compartment_size_distribution_by_subcomp(df, list(itertools.chain(*axs)))
    plt.tight_layout()

    outname = outprefix.parent / (outprefix.name + "_size_distr_by_subcomp")
    fig.savefig(outname.with_suffix(".svg"))
    fig.savefig(outname.with_suffix(".png"), dpi=300)
    plt.close()

    fig, axs = plt.subplots(1, 3, figsize=(4 * 6.4, 6.4))
    plot_compartment_size_distribution_by_condition(df, axs)
    outname = outprefix.parent / (outprefix.name + "_size_distr_by_condition")
    fig.savefig(outname.with_suffix(".svg"))
    fig.savefig(outname.with_suffix(".png"), dpi=300)
    plt.close()


if __name__ == "__main__":
    main()
