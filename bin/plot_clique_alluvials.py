#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import pathlib
import subprocess as sp
import sys
import tempfile
from typing import Union

import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bedgraph",
        type=pathlib.Path,
        nargs="+",
        help="Path to the TAD clique bedGraph (cliques should be called using same TAD annotation).",
    )

    cli.add_argument(
        "--path-to-plotting-script",
        type=pathlib.Path,
        default="plot_clique_alluvials.r",
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
        help="Skip plot generation. Only output a TSV with the data after aggregation.",
    )

    cli.add_argument("--labels", required=True, help="Comma-separated list of condition labels.")

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


def import_data(path_to_dfs: list, labels: list) -> pd.DataFrame:
    dfs = []

    for path, label in zip(path_to_dfs, labels):
        dfs.append(pd.read_table(path).drop(columns="tad_ids").set_index("name").rename(columns={"size": label}))

    return functools.reduce(
        lambda left, right: pd.merge(left, right, how="outer", left_index=True, right_index=True),
        dfs,
    )


def group_and_sort_clique_sizes(df: pd.DataFrame) -> pd.DataFrame:
    df = df.applymap(lambda x: 0 if x <= 2 else x)
    df = df.groupby(df.columns.tolist(), as_index=False).size()

    cols = [col for col in df.columns if col != "size"]
    return df.sort_values(
        by=cols,
        ignore_index=True,
        kind="stable",
    ).astype(int)


def plot_alluvial_with_r(
    path_to_input: pathlib.Path,
    output_name: pathlib.Path,
    path_to_rscript: pathlib.Path,
    label_to_highlight: Union[str, None] = None,
) -> None:
    cmd = [
        str(path_to_rscript),
        "--cliques",
        str(path_to_input),
        "--outprefix",
        str(output_name.with_suffix("")),
    ]
    if label_to_highlight is not None:
        cmd.append(f"--highlight_label={label_to_highlight}")

    try:
        sp.run(cmd, capture_output=True, text=True).check_returncode()
    except sp.CalledProcessError as e:
        cmd = " ".join((x for x in cmd))
        print(
            f'Process "{cmd}" terminated with exit code {e.returncode}.',
            file=sys.stderr,
        )
        if len(e.stdout) != 0:
            print(f"STDOUT:\n{e.stdout}", file=sys.stderr)
        if len(e.stderr) != 0:
            print(f"STDERR:\n{e.stderr}", file=sys.stderr)

        raise e


def make_alluvial_plot(
    path_to_tsv: pathlib.Path,
    label_to_highlight: str,
    outname: pathlib.Path,
    path_to_plotting_script: pathlib.Path,
    overwrite_existing: bool,
) -> None:
    if outname.exists() and not overwrite_existing:
        raise RuntimeError(f"Refusing to overwrite file {outname}. Pass --force to overwrite existing file(s).")

    plot_alluvial_with_r(path_to_tsv, outname, path_to_plotting_script, label_to_highlight)


def get_labels_from_df(df: pd.DataFrame) -> list:
    return pd.concat([df[col] for col in df.columns if col != "size"]).astype(int).unique().tolist()


def main():
    args = vars(make_cli().parse_args())
    df = group_and_sort_clique_sizes(import_data(args["bedgraph"], args["labels"].split(",")))

    output_prefix = pathlib.Path(args["output_prefix"])
    output_prefix.parent.mkdir(exist_ok=True)

    df.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)
    if args["no_plotting"]:
        return

    with tempfile.NamedTemporaryFile() as tmpdata:
        df.to_csv(tmpdata.name, sep="\t", index=False)

        for label in get_labels_from_df(df):
            outname = output_prefix.parent / f"{output_prefix.name}_{label}.svg"
            make_alluvial_plot(
                pathlib.Path(tmpdata.name),
                label,
                outname,
                args["path_to_plotting_script"],
                args["force"],
            )


if __name__ == "__main__":
    main()
