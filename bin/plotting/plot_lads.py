#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import gzip
import pathlib
import re
import subprocess as sp
import tempfile
from typing import List, Tuple

import bioframe as bf
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "--lads",
        required=True,
        nargs="+",
        type=pathlib.Path,
        help="Path to LADs in BED format.",
    )
    cli.add_argument(
        "--subcomps-pca",
        required=True,
        type=pathlib.Path,
        help="Path to the subcompartment bedgraph produced by dcHiC.",
    )

    cli.add_argument("--output", required=True, type=pathlib.Path, help="Path to output")

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    cli.add_argument(
        "--region",
        type=str,
        default="chr2:100000000-150000000",
        help="Region to use for plotting (UCSC format).",
    )

    return cli


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def init_ini(ini: pathlib.Path):
    args = ["[x-axis]", "where = top", "[spacer]", "height = 0.05"]

    with open(ini, "a") as f:
        print("\n".join(args), file=f)


def add_lads_to_ini(ini: pathlib.Path, *lads: pathlib.Path):
    for i, lad in enumerate(lads):
        args = [
            f"[lads_{i:03d}]",
            f"file = {lad}",
            "file_type = bed",
            "display = collapsed",
            "labels = false",
            "[spacer]",
            "height = 0.05",
        ]

        with open(ini, "a") as f:
            print("\n".join(args), file=f)


def add_compartments_to_ini(ini: pathlib.Path, *compartments: pathlib.Path, lb: float, ub: float):
    assert lb <= ub
    for i, comp in enumerate(compartments):
        args = [
            f"[compartments_{i:03d}]",
            f"file = {comp}",
            "file_type = bedgraph",
            f"min_value = {lb}",
            f"max_value = {ub}",
            "color = red",
            "negative_color = blue",
            "height = 1.0",
            "[spacer]",
            "height = 0.05",
        ]

        with open(ini, "a") as f:
            print("\n".join(args), file=f)


def import_pca(pca: pathlib.Path) -> List[pd.DataFrame]:
    df = (
        pd.read_table(pca)
        .rename(columns={"chr": "chrom"})
        .drop(columns=["replicate_wt", "sample_maha", "pval", "padj", "dist_clust"])
        .set_index(["chrom", "start", "end"])
    )

    return [df[col].to_frame() for col in df]


def import_lads(lads: List[pathlib.Path]) -> List[pd.DataFrame]:
    return [bf.read_table(doms, schema="bed") for doms in lads]


def compute_bedgraph_bounds(dfs: List[pd.DataFrame], region: str) -> Tuple[float, float]:
    lb = np.inf
    ub = -np.inf

    for df in dfs:
        col = df.columns[-1]
        lb = min(lb, bf.select(df.reset_index(), region)[col].dropna().min())
        ub = max(ub, bf.select(df.reset_index(), region)[col].dropna().max())

    return lb, ub


def extract_attribute_gtf(data: pd.Series, key: str) -> List[str]:
    pattern = re.compile(rf"{key} \"(.*?)\";")

    return data.str.extract(pattern)


def filter_gtf(input: pathlib.Path, output: pathlib.Path):
    with gzip.open(input) as fi, open(output, "w") as fo:
        for line in fi:
            line = line.decode("utf-8")
            if "#" in line:
                continue
            if "protein_coding" in line:
                print(line, file=fo)


def extract_gene_types(df: pd.DataFrame, gene_type_key="gene_type"):
    return np.sort(df[gene_type_key].unique())


def run_pygenometracks(ini: pathlib.Path, region: str, output: pathlib.Path):
    args = [
        "pyGenomeTracks",
        "--tracks",
        str(ini),
        "--region",
        region,
        "--outFileName",
        str(output),
    ]
    sp.check_call(args)


def main():
    args = vars(make_cli().parse_args())

    if not args["force"]:
        handle_path_collisions(args["output"])

    lads = import_lads(args["lads"])
    pcas = import_pca(args["subcomps_pca"])

    compartments_lb, compartments_ub = compute_bedgraph_bounds(pcas, args["region"])
    with tempfile.TemporaryDirectory() as tmpdir:
        ini = pathlib.Path(tmpdir) / "spec.ini"
        init_ini(ini)

        for i, df in enumerate(lads):
            path = pathlib.Path(tmpdir) / f"lads_{i:03d}.bed"
            df.to_csv(path, sep="\t", index=False, header=False)
            add_lads_to_ini(ini, path)

        for i, df in enumerate(pcas):
            path = pathlib.Path(tmpdir) / f"pcas_{i:03d}.bedgraph"
            df.to_csv(path, sep="\t", index=True, header=False)
            add_compartments_to_ini(ini, path, lb=compartments_lb, ub=compartments_ub)

        run_pygenometracks(ini, args["region"], args["output"])


if __name__ == "__main__":
    main()
