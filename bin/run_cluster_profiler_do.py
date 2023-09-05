#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import pathlib
from typing import Any

import matplotlib.pyplot as plt
import numpy.typing as npt
import pandas as pd
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def probability(s) -> float:
        if 0.0 <= (x := float(s)) <= 1:
            return x

        raise RuntimeError("Not a valid probability")

    cli.add_argument(
        "deg-table",
        type=pathlib.Path,
        help="Path to the table of differentially expressed genes.",
    )
    cli.add_argument(
        "output-prefix",
        type=pathlib.Path,
        help="Prefix path to use for output.",
    )

    cli.add_argument("--type", type=str, default="ora", choices={"ora", "gse"})

    cli.add_argument(
        "--pvalue",
        type=probability,
        default=0.01,
        help="P-value cutoff used to determine differentially expressed genes.",
    )

    cli.add_argument(
        "--qvalue-cp",
        type=probability,
        default=0.05,
        help="Q-value cutoff passed to clusterProfiler.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def import_de_genes(path_to_tsv):
    df = pd.read_table(path_to_tsv)
    gene_ids = df["id"].str.replace(r"\.\d+$", "", regex=True)

    with (ro.default_converter + pandas2ri.converter).context():
        df["entrez_id"] = annotation_dbi.mapIds(
            org_hs_eg_db.org_Hs_eg_db, keys=gene_ids, keytype="ENSEMBL", column="ENTREZID", multiVals="first"
        )

        # Replace R nans
        df["entrez_id"] = df["entrez_id"].apply(
            lambda val: pd.NA if isinstance(val, rpy2.rinterface_lib.sexp.NACharacterType) else val
        )

    return df.dropna()


def collect_de_gene_ids(df: pd.DataFrame, pvalue: float) -> npt.NDArray[str]:
    assert 0.0 <= pvalue <= 1.0
    de_genes = df.loc[df["svalue"] < pvalue, "entrez_id"]

    return de_genes.astype(str).to_numpy()


def prepare_for_gse(df: pd.DataFrame, pvalue: float) -> pd.DataFrame:
    assert 0.0 <= pvalue <= 1.0
    df = df[df["svalue"] < pvalue]
    df["log2FoldChange"] = df["log2FoldChange"].abs()
    df = df.sort_values("log2FoldChange", ascending=False)

    return df[["entrez_id", "log2FoldChange"]].rename(columns={"entrez_id": "id"})


def enrich_do_terms(de_genes: npt.NDArray[str], qval: float) -> Any:
    return dose.enrichDO(
        gene=rpy2.robjects.StrVector(de_genes),
        ont="DO",
        pAdjustMethod="BH",
        qvalueCutoff=qval,
        readable=False,
    )


def gse_do_terms(genes: pd.DataFrame, qval: float):
    genes_r = rpy2.robjects.FloatVector(genes["log2FoldChange"])
    genes_r.names = genes["id"].tolist()
    return dose.gseDO(gene=genes_r, minGSSize=10, pAdjustMethod="BH", pvalueCutoff=qval, nPermSimple=100000)


def make_dotplot(edo, output_name: pathlib.Path, showCategory: int = 15):
    try:
        enrichplot.dotplot(edo, showCategory=showCategory)
        ggplot.ggsave(str(output_name))
    except rpy2.rinterface_lib.embedded.RRuntimeError:
        # Save an empty plot
        fig, ax = plt.subplots(1, 1)
        fig.savefig(output_name)


def main():
    args = vars(make_cli().parse_args())

    if not args["force"]:
        handle_path_collisions(args["output-prefix"].with_suffix(".rds"), args["output-prefix"].with_suffix(".pdf"))

    df = import_de_genes(args["deg-table"])

    analysis_type = args["type"]
    if analysis_type == "ora":
        de = collect_de_gene_ids(df, args["pvalue"])
        robj = enrich_do_terms(de, args["qvalue_cp"])
    else:
        df_gse = prepare_for_gse(df, args["pvalue"])
        robj = gse_do_terms(df_gse, args["qvalue_cp"])

    ro.r["saveRDS"](robj, str(args["output-prefix"].with_suffix(".rds")), compress="xz")
    make_dotplot(robj, args["output-prefix"].with_suffix(".pdf"))


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)

    def log_consoleprint(msg):
        msg = msg.strip()
        if len(msg) != 0:
            logging.info("[R] %s", msg)

    def log_warnerror(msg):
        msg = msg.strip()
        if len(msg) != 0:
            logging.warning("[R] %s", msg)

    rpy2.rinterface_lib.callbacks.consolewrite_print = log_consoleprint
    rpy2.rinterface_lib.callbacks.consolewrite_warnerror = log_warnerror


if __name__ == "__main__":
    setup_logger()

    annotation_dbi = importr("AnnotationDbi")
    clusterProfiler = importr("clusterProfiler")
    dose = importr("DOSE")
    enrichplot = importr("enrichplot")
    ggplot = importr("ggplot2")
    org_hs_eg_db = importr("org.Hs.eg.db")
    main()
