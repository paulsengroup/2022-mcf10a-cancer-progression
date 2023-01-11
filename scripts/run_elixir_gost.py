#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import sys
import typing
from typing import List, Tuple, Union

import pandas as pd
from gprofiler import GProfiler


def make_cli():
    cli = argparse.ArgumentParser(
        "Read a TSV with the list of differentially expressed genes, then perform functional enrichment analysis using g:GOSt https://biit.cs.ut.ee/gprofiler/gost."
    )

    def valid_probability(arg):
        n = float(arg)
        if 0.0 <= n <= 1.0:
            return n

        raise RuntimeError("Not a valid probability")

    cli.add_argument(
        "--id-column-label",
        type=str,
        default="id",
        help="Name of the column containing gene/transcript IDs.",
    )

    cli.add_argument(
        "--lfc-cutoffs",
        type=float,
        nargs=2,
        default=[-1.0, 1.0],
        help="Log-fold change cutoff (lower and upper-bound).",
    )

    cli.add_argument(
        "--pval-cutoff",
        type=valid_probability,
        default=0.05,
        help="Adjusted pvalue cut-off.",
    )

    cli.add_argument(
        "--use-nonde-as-background",
        action="store_true",
        default=False,
        help="Use non differentially expressed genes as background.",
    )

    cli.add_argument(
        "--gprofiler-api-url",
        type=str,
        default="https://biit.cs.ut.ee/gprofiler_archive3/e105_eg52_p16",
    )

    return cli


def import_data(column_id: str, lfc_cutoffs: Tuple[float, float], pval_cutoff: float) -> pd.DataFrame:
    df = pd.read_table(sys.stdin).reset_index().rename(columns={"index": "id"}).dropna()
    df["significant"] = (~df["log2FoldChange"].between(*lfc_cutoffs)) & (df["padj"] <= pval_cutoff)

    return df[[column_id, "significant"]]


def run_gost(queries: List[str], background: Union[List[str], None], api_base_url: str, **kwargs) -> pd.DataFrame:
    if background is not None:
        kwargs["background"] = background

    kwargs["query"] = queries

    gp = GProfiler(return_dataframe=True, base_url=api_base_url)
    df = gp.profile(organism="hsapiens", **kwargs)
    df = typing.cast(pd.DataFrame, df)

    df["parents"] = df["parents"].apply(lambda x: ";".join(x))

    return df


def main():
    args = vars(make_cli().parse_args())

    df = import_data(args["id_column_label"], args["lfc_cutoffs"], args["pval_cutoff"])

    queries = df.loc[df["significant"], "id"].tolist()
    if args["use_nonde_as_background"]:
        background = df.loc[~df["significant"], "id"].tolist()
    else:
        background = None

    df = run_gost(queries, background, args["gprofiler_api_url"])
    df.to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
