#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import pathlib
from collections import defaultdict
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
import rpy2
import rpy2.robjects as ro
from rpy2.robjects import Formula, pandas2ri
from rpy2.robjects.packages import importr


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def non_negative_float(s) -> float:
        if (x := float(s)) >= 0:
            return x

        raise RuntimeError("Not a positive float")

    cli.add_argument(
        "count-table",
        nargs="+",
        type=pathlib.Path,
        help="Path to gene count table produced by nfcore/rnaseq.",
    )
    cli.add_argument(
        "design-table",
        type=pathlib.Path,
        help="Path to a TSV with the experiment design.\n"
        "Expected columns:\n"
        " - sample_id\n"
        " - new_sample_id - New sample IDs. Leave blank to keep the original sample IDs\n"
        " - condition - Condition ID (see DESeq2 docs)\n"
        " - seq_type - Sequence type (see DESeq2 docs)\n"
        " - contrast - Condition(s) to use as contrast",
    )

    cli.add_argument("output-folder", type=pathlib.Path, help="Path to folder where output files will be stored.")

    cli.add_argument(
        "--lfc-thresh",
        type=non_negative_float,
        default=0.0,
        help="Log2FoldChange threshold (see DESeq2 docs for more details).",
    )

    cli.add_argument(
        "--lfc-shrinkage-method",
        type=str,
        default="apeglm",
        choices={"apeglm", "ashr", "normal"},
        help="Log2FoldChange shrinkage method (see DESeq2 docs for more details).",
    )

    cli.add_argument(
        "--min-counts",
        type=int,
        default=10,
        help="Discard genes with less than --min-counts counts across replicates.",
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


def import_count_table(path_to_tsvs: List[pathlib.Path], round_counts: bool = True) -> pd.DataFrame:
    df = None
    # index_col are usually [tx, gene_name] or [gene_id, gene_name].
    # Any name should work, as long as the first column contains some kind of unique ID
    for p in path_to_tsvs:
        if df is None:
            df = pd.read_table(p, index_col=[0, 1])
            continue

        df = df.merge(
            pd.read_table(p, index_col=[0, 1]),
            left_index=True,
            right_index=True,
            how="outer",
        )

    if round_counts:
        df = df.round().astype(int)

    return df.sort_index()


def import_design_table(path_to_design: pathlib.Path) -> pd.DataFrame:
    df = pd.read_table(path_to_design)
    assert set(df.columns.tolist()) == {
        "sample_id",
        "new_sample_id",
        "condition",
        "seq_type",
        "contrast",
    }

    mask = df["new_sample_id"].isna()
    df.loc[mask, "new_sample_id"] = df.loc[mask, "sample_id"]

    return df


def import_data(path_to_tsvs: List[pathlib.Path], path_to_design: pathlib.Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    design = import_design_table(path_to_design)
    counts = import_count_table(path_to_tsvs)

    counts = sort_and_rename_columns(counts, design)
    design["sample_id"] = design["new_sample_id"]
    design = design.drop(columns="new_sample_id").set_index("sample_id")

    return counts, design


def get_condition_to_control_mappings(design: pd.DataFrame) -> Dict[str, List[str]]:
    contrasts = design["contrast"].transform(lambda s: None if pd.isna(s) else s.split(","))

    mappings = defaultdict(set)

    for condition, contrasts in zip(design["condition"].tolist(), contrasts.tolist()):
        if contrasts is not None:
            mappings[condition].update(contrasts)

    return {k: list(v) for k, v in mappings.items()}


def sort_and_rename_columns(counts: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    mappings = {k: v for _, (k, v) in design[["sample_id", "new_sample_id"]].iterrows() if k != v}
    return counts[design["sample_id"].tolist()].rename(columns=mappings)


def filter_lowly_expressed_genes(counts: pd.DataFrame, design: pd.DataFrame, min_counts: int) -> pd.DataFrame:
    if min_counts <= 0:
        return counts

    drop = np.full(len(counts), True, dtype=bool)
    for condition, samples in design.groupby("condition"):
        df = counts[samples.index.tolist()]
        drop &= df.sum(axis="columns") < min_counts

    logging.info("Dropping %d/%d genes with very low read count", np.sum(drop), len(drop))
    return counts[~drop]


def generate_deseq_dds(counts: pd.DataFrame, design: pd.DataFrame) -> Any:
    with (ro.default_converter + pandas2ri.converter).context():
        counts_rdf = ro.conversion.get_conversion().py2rpy(counts.droplevel(level=1, axis="index"))
        design_rdf = ro.conversion.get_conversion().py2rpy(design[["condition", "seq_type"]])
        assert list(ro.r["colnames"](counts_rdf)) == list(ro.r["row.names"](design_rdf))

        return deseq2.DESeqDataSetFromMatrix(countData=counts_rdf, colData=design_rdf, design=Formula("~condition"))


def deseq2_get_results(
    dds_r,
    contrast: str,
    condition: str,
    lfc_thresh: float,
    lfc_shrinkage_method: str,
) -> Tuple[pd.DataFrame, Any]:
    with (ro.default_converter + pandas2ri.converter).context():
        results_r = deseq2.lfcShrink(
            dds_r,
            coef=f"condition_{condition}_vs_{contrast}",
            lfcThreshold=lfc_thresh,
            type=lfc_shrinkage_method,
            parallel=False,
        )

        df = ro.conversion.get_conversion().rpy2py(ro.r["as.data.frame"](results_r))

    # Deal with case where lfc_thresh = 0
    if "svalue" not in df:
        df["svalue"] = df["padj"]

    # See https://github.com/stephens999/ashr/issues/123
    # and https://github.com/stephens999/ashr/issues/117
    df["svalue"] = np.maximum(df["svalue"], 0)
    df.insert(0, "contrast", contrast)
    df.insert(1, "condition", condition)
    df.index.name = "id"

    return df, results_r


def relevel_dds(dds: Any, column: str, ref: str) -> Any:
    relevel_dds_r = ro.r(
        f"""
        function(dds, column, ref) {{
            dds${column} <- relevel(dds${column}, ref=ref)
            return(dds)
        }}
        """
    )
    return relevel_dds_r(dds, column, ref)


def main():
    args = vars(make_cli().parse_args())

    outdir = args["output-folder"]

    counts, design = import_data(args["count-table"], args["design-table"])
    counts = filter_lowly_expressed_genes(counts, design, args["min_counts"])

    assert len(counts.index.names) == 2
    index1, index2 = counts.index.names

    condition_to_control_mappings = get_condition_to_control_mappings(design)
    assert len(condition_to_control_mappings) != 0

    dds_r = generate_deseq_dds(counts, design)

    outdir.mkdir(parents=True, exist_ok=True)
    for condition, contrasts in condition_to_control_mappings.items():
        for contrast in contrasts:
            assert condition != contrast
            logging.info("Processing %s vs %s...", contrast, condition)

            output_prefix = outdir / f"{contrast}_vs_{condition}"
            output_tsv = output_prefix.with_suffix(".tsv.gz")
            output_rds = output_prefix.with_suffix(".res.rds")
            if not args["force"]:
                handle_path_collisions(output_tsv, output_rds)

            dds_r = relevel_dds(dds_r, column="condition", ref=contrast)
            dds_r = deseq2.DESeq(dds_r, parallel=False)

            results, results_r = deseq2_get_results(
                dds_r,
                contrast,
                condition,
                lfc_thresh=args["lfc_thresh"],
                lfc_shrinkage_method=args["lfc_shrinkage_method"],
            )

            # Add back gene names to DF
            results = results.merge(
                counts.index.to_frame().set_index(index1),
                how="left",
                left_index=True,
                right_index=True,
            )

            results.insert(0, index2, results.pop(index2))

            results.to_csv(output_tsv, sep="\t", index=True)
            ro.r["saveRDS"](results_r, str(output_rds), compress="xz")

    output_rds = outdir / "dds.rds"
    if not args["force"]:
        handle_path_collisions(output_rds)
    ro.r["saveRDS"](dds_r, str(output_rds), compress="xz")

    output_sessioninfo = outdir / "r_sessioninfo.txt"
    if not args["force"]:
        handle_path_collisions(output_sessioninfo)

    with open(output_sessioninfo, "w") as f:
        print(utils.sessionInfo(), file=f)


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

    base = importr("base")
    utils = importr("utils")
    stats = importr("stats")
    deseq2 = importr("DESeq2")

    main()
