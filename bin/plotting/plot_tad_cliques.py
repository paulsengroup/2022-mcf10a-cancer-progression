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
import matplotlib as mpl
import numpy as np
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "--cliques",
        required=True,
        nargs="+",
        type=pathlib.Path,
        help="Path to cliques produced by robomics/call_tad_cliques.",
    )
    cli.add_argument(
        "--clique-doms",
        required=True,
        nargs="+",
        type=pathlib.Path,
        help="Path to domains produced by robomics/call_tad_cliques.",
    )
    cli.add_argument(
        "--domains",
        required=True,
        type=pathlib.Path,
        help="Path to a BED file with TADs used as input when running robomics/call_tad_cliques.",
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


def add_genes_to_ini(ini: pathlib.Path, gtf: pathlib.Path):
    args = [
        "[genes]",
        f"file = {gtf}",
        "height = 3",
        "file_type = gtf",
        "merge_transcripts = true",
        "merge_overlapping_exons = true",
        "labels = false",
        "labels_in_margin = true",
        "preferred_name = gene_name",
        "[spacer]",
        "height = 0.05",
    ]

    with open(ini, "a") as f:
        print("\n".join(args), file=f)


def add_tads_to_ini(ini: pathlib.Path, *tads: pathlib.Path):
    for i, tad in enumerate(tads):
        args = [
            f"[tads_{i:03d}]",
            f"file = {tad}",
            "file_type = bed",
            "display = collapsed",
            "color = bed_rgb",
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


def map_max_clique_size_to_domains(cliques: pd.DataFrame, domains: pd.DataFrame) -> pd.DataFrame:
    df = cliques.copy()

    df["tad_ids"] = df["tad_ids"].str.split(",")
    df["size"] = df.apply(lambda x: len(x["tad_ids"]), axis="columns")
    df = df.explode("tad_ids")
    ids = df.groupby("tad_ids")[["size"]].max()
    ids.index = ids.index.astype(int)

    df = domains.merge(ids, left_on="name", right_index=True).fillna(1)
    df["score"] = df["size"].astype(int)
    df["name"] = df["score"].astype(str)
    df["strand"] = "."

    return df[bf.SCHEMAS["bed6"]]


def color_domains(domains: pd.DataFrame, color1="forestgreen", color2="royalblue") -> pd.DataFrame:
    df = domains.copy()
    df["thickStart"] = df["start"]
    df["thickEnd"] = df["end"]
    df["rgb"] = ",".join(f"{round(n * 255):.0f}" for n in mpl.colors.to_rgb(color1))
    idx = np.arange(len(df)) % 2 == 0
    df.loc[idx, "rgb"] = ",".join(f"{round(n * 255):.0f}" for n in mpl.colors.to_rgb(color2))

    return df[bf.SCHEMAS["bed9"]]


def import_pca(pca: pathlib.Path) -> List[pd.DataFrame]:
    df = (
        pd.read_table(pca)
        .rename(columns={"chr": "chrom"})
        .drop(columns=["replicate_wt", "sample_maha", "pval", "padj", "dist_clust"])
        .set_index(["chrom", "start", "end"])
    )

    return [df[col].to_frame() for col in df]


def import_tads(
    clique_domains: List[pathlib.Path],
    cliques: List[pathlib.Path],
    domains: pathlib.Path,
) -> List[pd.DataFrame]:
    tads = []
    df1 = bf.read_table(domains, schema="bed3")
    df1["name"] = 1
    df1["score"] = 1
    for doms, cliques in zip(clique_domains, cliques):
        df2 = pd.read_table(cliques)
        df3 = bf.read_table(doms, schema="bed")

        df = map_max_clique_size_to_domains(df2, df3)
        df = pd.concat([df, df1]).drop_duplicates(subset=["chrom", "start", "end"], keep="first").reset_index()

        df = color_domains(bf.sort_bedframe(df))
        tads.append(df)

    return tads


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

    assert len(args["cliques"]) == len(args["clique_doms"])

    tads = import_tads(args["clique_doms"], args["cliques"], args["domains"])
    pcas = import_pca(args["subcomps_pca"])

    compartments_lb, compartments_ub = compute_bedgraph_bounds(pcas, args["region"])
    with tempfile.TemporaryDirectory() as tmpdir:
        ini = pathlib.Path(tmpdir) / "spec.ini"
        init_ini(ini)

        path = pathlib.Path(tmpdir) / "genes.gtf"
        filter_gtf(args["gtf"], path)
        add_genes_to_ini(ini, path)

        for i, df in enumerate(tads):
            path = pathlib.Path(tmpdir) / f"tads_{i:03d}.bed"
            df.to_csv(path, sep="\t", index=False, header=False)
            add_tads_to_ini(ini, path)

        for i, df in enumerate(pcas):
            path = pathlib.Path(tmpdir) / f"pcas_{i:03d}.bedgraph"
            df.to_csv(path, sep="\t", index=True, header=False)
            add_compartments_to_ini(ini, path, lb=compartments_lb, ub=compartments_ub)

        run_pygenometracks(ini, args["region"], args["output"])


if __name__ == "__main__":
    main()
