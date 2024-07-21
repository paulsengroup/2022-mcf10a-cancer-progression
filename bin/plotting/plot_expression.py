#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import gzip
import pathlib
import re
import shutil
import subprocess as sp
import tempfile
from typing import List, Tuple

import bioframe as bf
import matplotlib as mpl
import numpy as np
import pandas as pd
import pyBigWig


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "--deg-table",
        required=True,
        type=pathlib.Path,
        help="Path to table with the list of differentially expressed genes.",
    )
    cli.add_argument(
        "--subcomps-pca",
        required=True,
        type=pathlib.Path,
        help="Path to the subcompartment bedgraph produced by dcHiC.",
    )

    cli.add_argument("--gtf", required=True, type=pathlib.Path, help="Gene annotation in GTF format.")

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
        default="chr2:120000000-130000000",
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


def add_genes_to_ini(ini: pathlib.Path, bed: pathlib.Path):
    args = [
        "[genes]",
        f"file = {bed}",
        "height = 3",
        "labels = false",
        "color = bed_rgb",
        "file_type = bed",
        "merge_transcripts = true",
        "merge_overlapping_exons = true",
        "arrowhead_included = true",
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


def compute_bedgraph_bounds(dfs: List[pd.DataFrame], region: str) -> Tuple[float, float]:
    lb = np.inf
    ub = -np.inf

    for df in dfs:
        col = df.columns[-1]
        lb = min(lb, bf.select(df.reset_index(), region)[col].dropna().min())
        ub = max(ub, bf.select(df.reset_index(), region)[col].dropna().max())

    return lb, ub


def import_pca(pca: pathlib.Path) -> List[pd.DataFrame]:
    df = (
        pd.read_table(pca)
        .rename(columns={"chr": "chrom"})
        .drop(columns=["replicate_wt", "sample_maha", "pval", "padj", "dist_clust"])
        .set_index(["chrom", "start", "end"])
    )

    return [df[col].to_frame() for col in df]


def extract_attribute_gtf(data: pd.Series, key: str) -> List[str]:
    pattern = re.compile(rf"{key} \"(.*?)\";")

    return data.str.extract(pattern)


def filter_gtf(
    input: pathlib.Path,
    output: pathlib.Path,
    gene_types: Tuple[str] = ("protein_coding", "lncRNA"),
):
    if len(gene_types) == 0:
        shutil.copyfile(input, output)
        return

    if isinstance(gene_types, str):
        gene_types = [gene_types]

    pattern = re.compile("(" + ")|(".join(gene_types) + ")")
    with gzip.open(input) as fi, open(output, "w") as fo:
        for line in fi:
            line = line.decode("utf-8")
            if "#" in line:
                continue
            if pattern.search(line):
                print(line, file=fo)


def extract_gene_types(df: pd.DataFrame, gene_type_key="gene_type"):
    return np.sort(df[gene_type_key].unique())


def import_gtf(path_to_gtf: pathlib.Path) -> pd.DataFrame:
    df = bf.read_table(path_to_gtf, schema="gtf", comment="#")
    df = df[df["feature"] == "gene"]

    for key in ["gene_id", "gene_name", "gene_type"]:
        df[key] = extract_attribute_gtf(df["attributes"], key)

    return (
        df[(df["gene_type"] == "protein_coding") | (df["gene_type"] == "lncRNA")]
        .drop(columns="attributes")
        .set_index("gene_id")
        .sort_index()
    )


def classify_de_genes(
    path_to_deg: pathlib.Path,
    gtf: pd.DataFrame,
    lfc_cutoff: float = 0.5,
    pval_cutoff: float = 0.01,
) -> pd.DataFrame:
    df = pd.read_table(path_to_deg)[["id", "log2FoldChange", "svalue"]]
    df = (
        gtf.reset_index()
        .merge(df, left_on="gene_id", right_on="id", how="left")
        .dropna()
        .set_index(["gene_id", "gene_name", "gene_type"])
        .sort_index()
    )

    for col in ["start", "end"]:
        df[col] = df[col].astype(int)

    df["log2FoldChange"] = np.nan_to_num(df["log2FoldChange"])
    downreg = df[(df["log2FoldChange"] < -lfc_cutoff) & (df["svalue"] < pval_cutoff)].drop(
        columns=["log2FoldChange", "svalue"]
    )
    nonde = df[(df["log2FoldChange"].abs() < lfc_cutoff) | (df["svalue"] >= pval_cutoff)].drop(
        columns=["log2FoldChange", "svalue"]
    )
    upreg = df[(df["log2FoldChange"] > lfc_cutoff) & (df["svalue"] < pval_cutoff)].drop(
        columns=["log2FoldChange", "svalue"]
    )

    downreg["itemRgb"] = ",".join(f"{round(n * 255):.0f}" for n in mpl.colors.to_rgb("blue"))
    nonde["itemRgb"] = ",".join(f"{round(n * 255):.0f}" for n in mpl.colors.to_rgb("gray"))
    upreg["itemRgb"] = ",".join(f"{round(n * 255):.0f}" for n in mpl.colors.to_rgb("red"))

    df = pd.concat([downreg, nonde, upreg]).reset_index()

    df["thickStart"] = df["start"]
    df["thickEnd"] = df["end"]
    df["name"] = df["gene_name"]
    df["score"] = 0

    return df.reset_index()[bf.SCHEMAS["bed9"]]


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

    pcas = import_pca(args["subcomps_pca"])

    genes = classify_de_genes(args["deg_table"], import_gtf(args["gtf"]))

    compartments_lb, compartments_ub = compute_bedgraph_bounds(pcas, args["region"])
    with tempfile.TemporaryDirectory() as tmpdir:
        ini = pathlib.Path(tmpdir) / "spec.ini"
        init_ini(ini)

        path = pathlib.Path(tmpdir) / "genes.bed9"
        genes.to_csv(path, sep="\t", index=False, header=False)
        add_genes_to_ini(ini, path)

        for i, df in enumerate(pcas):
            path = pathlib.Path(tmpdir) / f"pcas_{i:03d}.bedgraph"
            df.to_csv(path, sep="\t", index=True, header=False)
            add_compartments_to_ini(ini, path, lb=compartments_lb, ub=compartments_ub)

        run_pygenometracks(ini, args["region"], args["output"])


if __name__ == "__main__":
    main()
