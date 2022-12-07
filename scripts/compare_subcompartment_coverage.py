#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import itertools
import math
import natsort
import bioframe as bf

import matplotlib.pyplot as plt
import pandas as pd

import seaborn as sns

import functools

import pathlib
from typing import List, Union


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
        "--genome-wide",
        action="store_true",
        default=False,
        help="Compute genome-wide subcompartment coverage.",
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


def compute_coverage(df: pd.DataFrame, filter: Union[str, None] = None) -> pd.DataFrame:

    subcompartments = get_compartment_labels_from_df(df)

    df = df.copy()
    if filter is not None:
        df = df[df["chrom"] == filter]

    df["span"] = df["end"] - df["start"]
    total_span = df["span"].sum()

    series = []
    for col in df.columns:
        if not col.endswith(".state"):
            continue

        cov = (
            (df[["span", col]].groupby(col).sum() / total_span)
            .rename(columns={"span": col})
            .T
        )
        series.append(cov)

    df = pd.concat(series, axis="index").fillna(0.0)

    for subcmp in subcompartments:
        if subcmp not in df.columns:
            df[subcmp] = 0.0

    df = df.T.sort_index(
        key=lambda x: pd.Index([get_compartment_ranks().get(y) for y in x]),
        kind="stable",
    )

    return df


def plot_coverage(
    df: pd.DataFrame, label: str, ax: plt.Axes, plot_legend: bool = False
) -> None:
    df.T.plot(kind="bar", stacked=True, ax=ax, legend=False, cmap="coolwarm")

    ax.set(
        title=f"Subcompartment coverage ({label})",
        xlabel="Conditions",
        ylabel="Relative coverage",
    )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
    if plot_legend:
        ax.legend(title=None)


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


def compute_coverage_genomewide(df: pd.DataFrame, outprefix: pathlib.Path, dpi: int = 300) -> None:
    coverage = compute_coverage(df)

    outname = outprefix.parent / (outprefix.name + "_genomewide")
    coverage.to_csv(outname.with_suffix(".tsv"), sep="\t")

    fig, ax = plt.subplots(1, 1)

    plot_coverage(coverage, "genome-wide", ax, plot_legend=True)
    plt.tight_layout()
    fig.savefig(outname.with_suffix(".svg"))
    fig.savefig(outname.with_suffix(".png"), dpi=dpi)
    plt.close()


def compute_coverage_by_chrom(df: pd.DataFrame, outprefix: pathlib.Path, dpi: int = 300) -> None:
    chroms = natsort.natsorted(df["chrom"].unique())
    plot_grid_size = int(math.ceil(math.sqrt(len(chroms))))
    fig, axs = plt.subplots(
        plot_grid_size,
        plot_grid_size,
        figsize=(6.4 * plot_grid_size, 6.4 * plot_grid_size),
    )

    dfs = []
    plot_legend = True
    for ax, chrom in zip(itertools.chain(*axs), chroms):
        coverage = compute_coverage(df, chrom)
        plot_coverage(coverage, chrom, ax, plot_legend=plot_legend)
        plot_legend = False
        coverage.insert(loc=0, column="chrom", value=pd.Series([chrom] * len(coverage)))
        dfs.append(coverage)

    outname = outprefix.parent / (outprefix.name + "_by_chrom")

    plt.tight_layout()
    fig.savefig(outname.with_suffix(".svg"))
    fig.savefig(outname.with_suffix(".png"), dpi=dpi)
    plt.close()

    pd.concat(dfs).to_csv(outname.with_suffix(".tsv"))


def main():
    args = vars(make_cli().parse_args())

    df = import_data(args["bedgraph"])

    if args["genome_wide"]:
        compute_coverage_genomewide(df, args["output-prefix"])
    else:
        compute_coverage_by_chrom(df, args["output-prefix"])


if __name__ == "__main__":
    main()
