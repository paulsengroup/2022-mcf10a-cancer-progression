#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import logging
import pathlib
import pickle
import sys
from typing import List, Tuple, Union

import hdbscan
import pandas as pd


def make_cli():
    def count_or_fraction(arg) -> float:
        n = float(arg)
        if n >= 1 and n.is_integer():
            return n
        if 0 <= n < 1:
            return n

        raise ValueError("Not a count or fraction")

    def comma_separated_string(s: str) -> Tuple[str]:
        return tuple(s.split(","))

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "tsv",
        nargs="+",
        type=pathlib.Path,
        help="TSV or BED file produced by annotate_domains_with_subcompartments.py.",
    )
    cli.add_argument(
        "--labels",
        type=comma_separated_string,
        help="Comma-separated list of labels used to label conditions. Should be in the same order as the input TSV(s).",
    )
    cli.add_argument(
        "--min-cluster-size",
        type=count_or_fraction,
        default=0.01,
        help="Minimum cluster size.\n"
        "Values in range [0, 1) are interpreted as fractions of the total number of datapoints.\n"
        "See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--min-samples",
        type=count_or_fraction,
        default=5,
        help="Minimum samples.\n"
        "Values in range [0, 1) are interpreted as fractions of the total number of datapoints.\n"
        "See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--cluster-selection-method",
        type=str,
        choices={"eom", "leaf"},
        default="leaf",
        help="See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--cluster-selection-epsilon",
        type=float,
        default=0.05,
        help="See HDBSCAN* documentation for more details.",
    )
    cli.add_argument(
        "--distance-metric",
        type=str,
        choices={"euclidean", "manhattan", "chebyshev"},
        default="euclidean",
        help="Distance metric used for clustering.",
    )
    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Path to use as output prefix (including parent folder(s) but without extension).",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


@functools.cache
def get_subcompartment_ranks() -> dict:
    compartment_labels = tuple(["B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3"])
    return {k: v for v, k in enumerate(compartment_labels)}


@functools.cache
def get_compartment_ranks() -> dict:
    compartment_labels = tuple(["B", "A"])
    return {k: v for v, k in enumerate(compartment_labels)}


def import_tsv(path_to_tsvs: List[pathlib.Path], labels: Union[List[str], None]) -> pd.DataFrame:
    assert labels is None or len(path_to_tsvs) == len(labels)
    stdin = {pathlib.Path("stdin"): sys.stdin, pathlib.Path("-"): sys.stdin}

    paths = ((stdin.get(path, path)) for path in path_to_tsvs)
    if labels is None:
        labels = [p.name.partition(".")[0] for p in path_to_tsvs]

    dfs = []
    for label, path in zip(labels, paths):
        df = pd.read_table(path)
        df["label"] = label
        dfs.append(df)

    return pd.concat(dfs)


def run_clustering(
    df: pd.DataFrame, dist_metric, min_cluster_size, min_samples, cluster_selection_method, cluster_selection_epsilon
) -> Tuple[pd.DataFrame, hdbscan.HDBSCAN]:
    cols = df.filter(regex=r"[AB]\d\.state").columns.tolist()
    m = df[cols].to_numpy()
    m = m / m.sum(axis=1)[:, None]

    clusterer = hdbscan.HDBSCAN(
        metric=dist_metric,
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method=cluster_selection_method,
        cluster_selection_epsilon=cluster_selection_epsilon,
    )
    clusterer.fit_predict(m)

    df = df.copy()
    df["cluster"] = clusterer.labels_
    df["cluster_pval"] = clusterer.probabilities_
    df["outlier_score"] = clusterer.outlier_scores_
    df[cols] = m

    return df, clusterer


def compute_min_cluster_size(df: pd.DataFrame, n: float) -> int:
    if n >= 1:
        return int(n)

    return int(round(len(df) * n))


def compute_min_samples(df: pd.DataFrame, n: float) -> int:
    return compute_min_cluster_size(df, n)


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def main():
    args = vars(make_cli().parse_args())

    cluster_selection_method = args["cluster_selection_method"]
    cluster_selection_epsilon = args["cluster_selection_epsilon"]

    out_prefix = args["output_prefix"]

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    df = import_tsv(args["tsv"], args.get("labels"))

    min_cluster_size = compute_min_cluster_size(df, args["min_cluster_size"])
    min_samples = compute_min_samples(df, args["min_samples"])

    df, clusterer = run_clustering(
        df,
        dist_metric=args["distance_metric"],
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        cluster_selection_method=cluster_selection_method,
        cluster_selection_epsilon=cluster_selection_epsilon,
    )

    outname = out_prefix.with_suffix(".clusters.tsv.gz")
    if not args["force"]:
        handle_path_collisions(outname)

    df.to_csv(outname, sep="\t", index=False, header=True)

    outname = out_prefix.with_suffix(".clusterer.pickle.xz")
    if not args["force"]:
        handle_path_collisions(outname)

    with open(outname, "wb") as f:
        pickle.dump(clusterer, f, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    setup_logger()
    main()
