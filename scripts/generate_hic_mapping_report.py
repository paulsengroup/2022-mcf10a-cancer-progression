#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import csv
import pathlib
import tarfile
import tempfile
import warnings
from typing import Tuple, Union

import matplotlib.pyplot as plt
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "reports",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more TAR files containing HiC mapping stats.",
    )
    cli.add_argument("output-prefix", type=pathlib.Path, help="Path to the output prefix.")
    cli.add_argument(
        "--sample-labels",
        type=str,
        help="Comma-separated list of sample labels to use instead of those inferred from file names.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    return cli


def generate_tar_member_list(tar: tarfile.TarFile, suffixes: Union[None, list] = None) -> Tuple[str, dict]:
    if suffixes is None:
        suffixes = ["allValidPairs.mergestat", "mpairstat", "mRSstat"]

    members = {}
    for member in tar.getmembers():
        for suffix in suffixes:
            if member.path.endswith(suffix):
                members[suffix] = member
                break

    if len(members) != len(suffixes):
        missing_files = ", ".join((s for s in suffixes if s not in members.keys()))
        raise RuntimeError(f"Unable to extract files with suffixes {missing_files} from TAR file {tar.name}")

    return tar.getmembers()[0].path, members


def import_data_from_tar(path: pathlib.Path) -> pd.DataFrame:
    dfs = []
    with tempfile.TemporaryDirectory() as tmpdir:
        with tarfile.open(path) as tar:
            name, members = generate_tar_member_list(tar)
            tar.extractall(tmpdir, members=members.values())

        for suffix, member in members.items():
            path = pathlib.Path(tmpdir) / member.path
            df = pd.read_table(path, header=None, index_col=0, usecols=range(2)).T
            assert len(df) == 1
            df["sample"] = name
            dfs.append(df.set_index("sample").T)

    df = pd.concat(dfs).T
    df.index.set_names("sample")
    return df


def merge_dfs(dfs: list, labels: Tuple[None, list] = None) -> pd.DataFrame:
    df = pd.concat(dfs)

    if labels is not None:
        assert len(df) == len(labels)
        df.index = pd.Index(labels)

    df.index.name = "sample"
    return df


def make_plot_contact_type(
    df: pd.DataFrame,
    relative: bool = False,
) -> plt.Figure:

    columns = {
        "cis_longRange": "cis > 20kbp",
        "cis_shortRange": "cis < 20kbp",
        "trans_interaction": "trans",
    }

    df1 = df[list(columns.keys())].rename(columns=columns)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        df1["filtered"] = df["valid_interaction"] - df["valid_interaction_rmdup"]

    if relative:
        df1 = df1.divide(df["valid_interaction"], axis="rows")
    else:
        df1 = df1 / 1.0e6

    table = False
    if not relative:
        table = df1.copy()
        table["total"] = df1.sum(axis="columns")
        table = table.round(2).T

    fig, ax = plt.subplots(1, 1)
    df1.plot(kind="bar", stacked=True, ax=ax, legend=False, table=table)

    if relative:
        title = "HiC relative interactions"
        ylabel = "Relative interactions"
        xlabel = "Samples"
        ax.set_xticklabels(ax.get_xticklabels(), rotation=30)
    else:
        title = "HiC interactions"
        ylabel = "Interactions (millions)"
        xlabel = ""
        ax.get_xaxis().set_ticks([])

    ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
    ax.legend(loc="lower right", title=None)

    return fig


def check_file_name_collisions(*paths: Union[str, pathlib.Path]) -> None:
    existing_files = []
    for path in paths:
        if path.exists():
            existing_files.append(path)

    if len(existing_files) != 0:
        existing_files = ", ".join((str(f) for f in existing_files))
        raise RuntimeError(
            f"Refusing to overwrite file(s) {existing_files}. Pass --force to overwrite existing file(s)."
        )


def main():
    args = vars(make_cli().parse_args())
    labels = args.get("sample_labels")

    if labels is not None:
        labels = labels.split(",")
        if (num_files := len(args["reports"])) != len(labels):
            raise RuntimeError(f"Expected {num_files} labels, found {len(labels)}")

    path_to_out_table = args["output-prefix"].with_suffix(".tsv")
    plot_out_prefix = args["output-prefix"]

    if not args["force"]:
        check_file_name_collisions(
            path_to_out_table,
            plot_out_prefix.with_suffix(".svg"),
            plot_out_prefix.with_suffix(".png"),
            plot_out_prefix.parent / (plot_out_prefix.name + "_relative.svg"),
            plot_out_prefix.parent / (plot_out_prefix.name + "_relative.png"),
        )

    df = merge_dfs([import_data_from_tar(path) for path in args["reports"]], labels)

    df["cis_trans_ratio"] = df["cis_interaction"] / df["trans_interaction"]
    df["cis_long_short_ratio"] = df["cis_longRange"] / df["cis_shortRange"]

    path_to_out_table.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path_to_out_table, sep="\t", quoting=csv.QUOTE_NONNUMERIC)

    fig = make_plot_contact_type(df)
    plt.tight_layout()
    fig.savefig(plot_out_prefix.with_suffix(".svg"))
    fig.savefig(plot_out_prefix.with_suffix(".png"), dpi=600)
    fig.clf()

    fig = make_plot_contact_type(df, relative=True)
    plt.tight_layout()
    fig.savefig(plot_out_prefix.parent / (plot_out_prefix.name + "_relative.svg"))
    fig.savefig(plot_out_prefix.parent / (plot_out_prefix.name + "_relative.png"), dpi=600)
    fig.clf()


if __name__ == "__main__":
    main()
