#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import collections
import functools
import itertools
import logging
import pathlib
import sys
import typing
from typing import Dict, Iterable, List, Tuple, Union

import bioframe as bf
import numpy as np
import numpy.typing as npt
import pandas as pd


def make_cli():
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "dchic-bedgraph",
        type=pathlib.Path,
        help="Path to the (sub)compartment bedGraph generated by dcHiC.",
    )
    cli.add_argument(
        "--domains",
        type=pathlib.Path,
        nargs="+",
        required=True,
        help="One or more domain files in BED format generated by robomics/call_tad_cliques.",
    )
    cli.add_argument(
        "--cliques",
        type=pathlib.Path,
        nargs="+",
        help="One or more clique file in TSV format generated by robomics/call_tad_cliques.",
    )

    cli.add_argument(
        "--domain-names",
        type=str,
        nargs="+",
        required=True,
        help="One or more condition name from the bedgraph generated by dcHiC.\n"
        "This is used to pair subcompartments labels from a given condition with cliques and domains.\n"
        "Thus, domain names should be listed in the same order as domain/clique files.",
    )

    cli.add_argument(
        "--mask",
        type=pathlib.Path,
        help="Path to a BED3 file with the genomic coordinates to be masked out.",
    )

    cli.add_argument(
        "-o",
        "--output-folder",
        type=pathlib.Path,
        required=True,
        help="Path to folder where to store output files.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    cli.add_argument(
        "--aggregate-subcompartments",
        action="store_true",
        default=False,
        help="Aggregate subcompartments to A/B compartments.",
    )

    return cli


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


@functools.cache
def get_subcompartment_ranks() -> dict:
    compartment_labels = tuple(["B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3"])
    return {k: v for v, k in enumerate(compartment_labels)}


@functools.cache
def get_compartment_ranks() -> dict:
    compartment_labels = tuple(["B", "A"])
    return {k: v for v, k in enumerate(compartment_labels)}


@functools.cache
def get_output_df_columns(aggregate_subcompartments: bool) -> List[str]:
    if aggregate_subcompartments:
        cols = [f"{s}.state" for s in get_compartment_ranks().keys()]
    else:
        cols = [f"{s}.state" for s in get_subcompartment_ranks().keys()]
    return cols + ["state.mode", "state.stderr", "state.mean_abs_dev"]


def states_contain_subcompartments(states: Iterable[str]) -> bool:
    unique_states = set(states)

    compartment_labels = set(list(get_compartment_ranks().keys()))
    return len(unique_states.intersection(compartment_labels)) == 0


def import_subcomps(
    path_to_bedgraph: pathlib.Path, labels: List[str], aggregate_subcompartments: bool
) -> List[pd.DataFrame]:
    df = pd.read_table(path_to_bedgraph)
    df1 = df[["chrom", "start", "end", "padj"]]

    dfs = []
    for label in labels:
        state_col = f"{label}.state"
        score_col = f"{label}.score"
        df2 = df1.copy()
        if aggregate_subcompartments:
            df2["state"] = df[state_col].str.extract(r"^([AB])\d+$")
            df2["score"] = df[score_col].map(lambda n: -1 if n < 0 else +1)
        else:
            df2["state"] = df[state_col]
            df2["score"] = df[score_col]

        dfs.append(df2)

    return dfs


def import_cliques(path_to_cliques: Union[pathlib.Path, None]) -> Union[pd.DataFrame, None]:
    if path_to_cliques is None:
        return None
    df = pd.read_table(path_to_cliques).rename(columns={"name": "clique"})
    assert df.columns.tolist() == ["clique", "tad_ids", "size"]
    df["tad_ids"] = df["tad_ids"].apply(lambda ids: [int(tid) for tid in ids.split(",")])
    return df.set_index("clique")


def compute_max_tad_clique_size(cliques: pd.DataFrame) -> pd.DataFrame:
    clique_sizes = {}
    clique_domains = {}

    for _, (tads, size) in cliques[["tad_ids", "size"]].iterrows():
        for tad in tads:
            if tad in clique_sizes:
                if size > clique_sizes[tad]:
                    clique_sizes[tad] = size
                    clique_domains[tad] = tads
            else:
                clique_sizes[tad] = size
                clique_domains[tad] = tads

    return pd.DataFrame({"tad": clique_sizes.keys(), "tad_ids": clique_domains.values(), "size": clique_sizes.values()})


def import_domains(path_to_domains: pathlib.Path) -> pd.DataFrame:
    return (
        pd.read_table(path_to_domains, names=bf.SCHEMAS["bed6"][0:4], usecols=list(range(0, 4)))
        .rename(columns={"name": "id"})
        .set_index("id")
    )


def overlap_domains_with_subcomps(subcomps: pd.DataFrame, domains: pd.DataFrame, mask: pd.DataFrame) -> pd.DataFrame:
    index_name = domains.index.name
    df = (
        bf.overlap(domains.reset_index(), subcomps, suffixes=("", "_"))
        .drop(columns=["chrom_", "start_", "end_"])
        .rename(columns=lambda c: c.rstrip("_"))
        .convert_dtypes()
        .set_index(index_name)
    )

    return bf.setdiff(df, mask)


def compute_state_distribution(df: pd.DataFrame) -> collections.Counter:
    return collections.Counter(df["state"].tolist())


def compute_modal_state(counts: collections.Counter) -> List[str]:
    mode = []
    mode_freq = 0
    for comp, count in counts.most_common():
        if count > mode_freq:
            mode = [comp]
            mode_freq = count
        elif count == mode_freq:
            mode.append(comp)
        else:
            break

    return mode


def compute_stderr(scores: npt.NDArray) -> float:
    return np.std(scores) / np.sqrt(len(scores))


def compute_mean_absolute_deviation(scores: npt.NDArray) -> float:
    return typing.cast(float, np.mean(np.abs(scores - np.mean(scores))))


def compute_state_stats(df: pd.DataFrame) -> Tuple:
    df_contains_subcompartments = states_contain_subcompartments(df["state"])
    ranks = get_subcompartment_ranks() if df_contains_subcompartments else get_compartment_ranks()

    if len(df) == 0:
        counts = np.zeros(len(ranks), dtype=int)
        return *counts, "NA", 0, 0

    counts = compute_state_distribution(df)
    mode = compute_modal_state(counts)

    counts = np.array([counts[subcmp] for subcmp in ranks.keys()])
    scores = df["state"].map(ranks).to_numpy()

    return (
        *counts,
        ",".join(mode),
        compute_stderr(scores),
        compute_mean_absolute_deviation(scores),
    )


def compute_stats_for_tad(key, grp):
    return list(key) + list(compute_state_stats(grp[["state", "score"]].dropna()))


def compute_state_stats_for_tads(df: pd.DataFrame) -> pd.DataFrame:
    data = []

    for key, grp in df.groupby(by=["chrom", "start", "end"]):
        data.append(compute_stats_for_tad(key, grp))

    df_has_subcompartments = states_contain_subcompartments(df["state"])
    cols = ["chrom", "start", "end"] + get_output_df_columns(aggregate_subcompartments=not df_has_subcompartments)
    return pd.DataFrame(data, columns=cols)


def compute_state_stats_for_clique(max_tad: int, domains: pd.DataFrame):
    cols = domains.filter(regex=r"[AB\d+]\.state$").columns.tolist()
    states = pd.DataFrame(data=[domains.loc[max_tad, cols]], columns=cols)
    states = pd.Series(
        itertools.chain.from_iterable(
            ([idx.removesuffix(".state")] * count for idx, count in states.sum(axis="index").items())
        ),
        name="state",
    )

    return [max_tad] + list(compute_state_stats(states.to_frame()))


def annotate_cliques(cliques: pd.DataFrame, domains: pd.DataFrame):
    data = []
    for _, (max_tad, tad_ids, size) in cliques.iterrows():
        data.append(compute_state_stats_for_clique(max_tad, domains) + [size])

    num_states = len(domains.filter(regex=r"^[AB\d]+\.state$").columns)

    return pd.DataFrame(
        data=data,
        columns=["clique"] + get_output_df_columns(aggregate_subcompartments=num_states <= 2) + ["clique_size"],
    )


def annotate_domains(subcomps: pd.DataFrame, domains: pd.DataFrame, mask: pd.DataFrame) -> pd.DataFrame:
    return compute_state_stats_for_tads(overlap_domains_with_subcomps(subcomps, domains, mask))


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())

    domain_files = args["domains"]
    clique_files = args["cliques"]
    if clique_files is None:
        clique_files = [None] * len(domain_files)

    aggregate_subcompartments = args["aggregate_subcompartments"]

    labels = args["domain_names"]
    if len(domain_files) != len(labels):
        raise RuntimeError(f"Expected {len(labels)} domain file(s), found {len(domain_files)}")

    if len(clique_files) != len(labels):
        raise RuntimeError(f"Expected {len(labels)} clique file(s), found {len(clique_files)}")

    if args["mask"] is not None:
        mask = pd.read_table(args["mask"], names=["chrom", "start", "end"])
    else:
        mask = pd.DataFrame([["a", 0, 0]], columns=["chrom", "start", "end"]).drop(0)

    subcomp_dfs = import_subcomps(args["dchic-bedgraph"], labels, aggregate_subcompartments)

    outfolder = args["output_folder"]
    outfolder.mkdir(parents=True, exist_ok=True)

    for label, subcomps, dom_bed, clique_tsv in zip(labels, subcomp_dfs, domain_files, clique_files):
        logging.info("[%s] annotating domains...", label)
        outfile = outfolder / dom_bed.name
        if not args["force"]:
            handle_path_collisions(outfile)

        doms = annotate_domains(subcomps, import_domains(dom_bed), mask)
        logging.info("[%s] writing domains to %s...", label, outfile)
        doms.to_csv(outfile, sep="\t", index=False, header=True, na_rep="nan")

        if clique_tsv is None:
            logging.info("[%s] clique file not available. SKIPPING clique annotation!", label)
            continue

        logging.info("[%s] annotating cliques...", label)
        outfile = outfolder / clique_tsv.name
        if not args["force"]:
            handle_path_collisions(outfile)

        cliques = import_cliques(clique_tsv)
        cliques = compute_max_tad_clique_size(cliques)
        cliques = annotate_cliques(cliques, doms)
        logging.info("[%s] writing cliques to %s...", label, outfile)
        cliques.to_csv(outfile, sep="\t", index=False, header=True, na_rep="nan")


if __name__ == "__main__":
    setup_logger()
    main()
