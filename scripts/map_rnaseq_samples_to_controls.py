#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import pathlib

import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "gene-count-table",
        type=pathlib.Path,
        help="TSV with the gene count table produced by nf-core/rnaseq.",
    )
    cli.add_argument(
        "sample-mappings",
        type=pathlib.Path,
        help="Path to a 2-column TSV with the control to sample mappings.",
    )
    cli.add_argument("output-dir", type=pathlib.Path, help="Path to output folder.")

    cli.add_argument(
        "--round",
        action="store_true",
        default=False,
        help="Round counts and output them as integers.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    return cli


def generate_sample_mappings(path_to_mappings: pathlib.Path) -> dict:
    df = pd.read_table(path_to_mappings)
    mappings = {}
    for control, samples in df.groupby("control"):
        mappings[control] = [control] + samples["sample"].tolist()

    return mappings


def import_counts(path_to_counts: pathlib.Path, round_counts: bool) -> pd.DataFrame:
    count_df = pd.read_table(path_to_counts, index_col=0)
    count_df.drop(columns=count_df.columns[0])
    count_df.index.name = "id"
    if round_counts:
        return count_df.round().apply(
            pd.to_numeric, errors="ignore", downcast="integer"
        )

    return count_df


def main():
    args = vars(make_cli().parse_args())

    path_to_gene_count_table = pathlib.Path(args["gene-count-table"])
    counts = import_counts(path_to_gene_count_table, args["round"])
    mappings = generate_sample_mappings(args["sample-mappings"])

    outdir = args["output-dir"]
    outdir.mkdir(exist_ok=True)

    for (i, (control, samples)) in enumerate(mappings.items()):
        pattern = "(" + ")|(".join((str(x) for x in samples)) + ")"
        outpath = outdir / pathlib.Path(
            f"{i:03d}_{path_to_gene_count_table.stem}_{control}.tsv"
        )

        if outpath.exists() and not args["force"]:
            raise RuntimeError(
                f"Refusing to overwrite file {outpath}. Pass --force to overwrite existing file(s)."
            )

        counts.filter(regex=pattern).to_csv(outpath, sep="\t")
        with open(f"{outpath}.contrast", "w") as f:
            print(control, file=f)


if __name__ == "__main__":
    main()
