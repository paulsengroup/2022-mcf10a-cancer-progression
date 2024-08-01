#!/usr/bin/env python3

# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import itertools
import logging
import pathlib
import re
import warnings
from collections import Counter
from typing import Iterable, List, Tuple, Union

import bioframe as bf
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import scipy.stats as ss
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm


def make_cli():
    cli = argparse.ArgumentParser()

    def existing_file(arg):
        if (path := pathlib.Path(arg)).exists():
            return path

        raise FileNotFoundError(f'File "{arg}" does not exists')

    def non_negative_float(s) -> float:
        if (x := float(s)) >= 0:
            return x

        raise RuntimeError("Not a non-negative float")

    def probability(s) -> float:
        if 0.0 <= (x := float(s)) <= 1.0:
            return x

        raise RuntimeError("Not a valid probability")

    def comma_separated_string(s: str) -> Tuple[str]:
        return tuple(s.split(","))

    cli.add_argument(
        "bedgraph",
        type=existing_file,
        help="Path to the (sub)compartment bedGraph generated by dcHiC.",
    )

    cli.add_argument(
        "de-table",
        type=existing_file,
        help="Path to a TSV with DE genes produced by run_deseq2.py.",
    )

    cli.add_argument(
        "gtf",
        type=existing_file,
        help="Path to a GTF with coordinates for genes found in the expression table.",
    )

    cli.add_argument(
        "output-prefix",
        type=pathlib.Path,
        help="Path to output prefix (including parent folder(s) but without extension).",
    )

    cli.add_argument(
        "--plot-type",
        type=str,
        choices={"boxplot", "pvalue-heatmap", "heatmap", "correlation"},
        default="correlation",
        help="Plot type.",
    )

    cli.add_argument(
        "--gene-types",
        type=comma_separated_string,
        default=["protein_coding", "lncRNA"],
        help="Comma-separated list of gene biotype(s).",
    )

    cli.add_argument(
        "--contrast",
        type=str,
        required=True,
        help="Name of the condition to use as contrast.",
    )

    cli.add_argument(
        "--condition",
        type=str,
        required=True,
        help="Name of the condition to use as treatment.",
    )

    cli.add_argument(
        "--lfc",
        type=non_negative_float,
        default=0.0,
        help="LFC cutoff. Ignored when --plot-type=correlation.",
    )

    cli.add_argument(
        "--padj",
        type=probability,
        default=0.01,
        help="LFC cutoff. Ignored when --plot-type=correlation.",
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


def handle_path_collisions(*paths: Union[str, pathlib.Path]) -> None:
    collisions = [p for p in paths if pathlib.Path(p).exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def import_subcomps(path_to_bedgraph: pathlib.Path, cond1: str, cond2: str) -> pd.DataFrame:
    """
    Read a bedgraph into a df with the following columns:
    chrom, start, end, padj, *.state
    where *.state stores the subcompartment label for the given interval.
    """
    logging.info("Reading subcompartments from %s...", path_to_bedgraph)
    df = pd.read_table(path_to_bedgraph).rename(columns={"chr": "chrom"}).set_index(["chrom", "start", "end"])
    logging.info("Imported %d intervals", len(df))
    for cond in (cond1, cond2):
        if f"{cond}.state" not in df:
            raise RuntimeError(f"{cond} not found in {path_to_bedgraph}")
    return df.filter(regex=".state$").rename(columns=lambda c: c.removesuffix(".state"))[[cond1, cond2]]


def extract_attribute_gtf(data: pd.Series, key: str) -> List[str]:
    pattern = re.compile(rf"{key} \"(.*?)\";")

    return data.str.extract(pattern)


def import_gtf(path_to_gtf: pathlib.Path) -> pd.DataFrame:
    logging.info("Importing genes from GTF...")
    df = bf.read_table(path_to_gtf, schema="gtf", comment="#")
    df = df[df["feature"] == "gene"]

    for key in ["gene_id", "gene_name", "gene_type"]:
        df[key] = extract_attribute_gtf(df["attributes"], key)

    return df.drop(columns="attributes").set_index("gene_id").sort_index()


def extract_gene_types(dfs: Iterable[pd.DataFrame], gene_type_key="gene_type"):
    return np.sort(np.unique(list(itertools.chain.from_iterable((df[gene_type_key].unique() for df in dfs)))))


def import_de_gene_table(path_to_tsv: pathlib.Path, contrast: str, condition: str) -> pd.DataFrame:
    df = pd.read_table(path_to_tsv).set_index("id")
    df.index.name = "gene_id"
    assert df["contrast"].isin([contrast]).any()
    assert df["condition"].isin([condition]).any()

    return df[(df["contrast"] == contrast) & (df["condition"] == condition)].drop(columns=["contrast", "condition"])


def get_modal_subcompartment(counter: Counter) -> Union[str, None]:
    """
    Find the most frequent subcompartment label.
    Returns None in case of ties.
    """
    if len(counter) == 1:
        return counter.most_common(1)[0][0]

    (s1, n1), (_, n2) = counter.most_common(2)
    if n1 == n2:
        return None
    return s1


def overlap_sucompartments_with_genes(subcomps: pd.DataFrame, de_table: pd.DataFrame, gtf: pd.DataFrame):
    labels = subcomps.columns.tolist()
    df = (
        bf.overlap(subcomps.reset_index(), gtf.reset_index(), how="inner", suffixes=("_", ""))
        .drop(columns=["chrom_", "start_", "end_"])
        .rename(columns=lambda c: c.rstrip("_"))
    )
    df = df[["gene_id", "gene_type"] + labels].groupby("gene_id").aggregate(Counter)

    df = df.applymap(get_modal_subcompartment)
    return de_table.merge(df, left_index=True, right_index=True)


def overlap_sucompartments_with_tss(subcomps: pd.DataFrame, de_table: pd.DataFrame, gtf: pd.DataFrame):
    gtf1 = gtf.copy().reset_index()
    gtf1["end"] = gtf1["start"] + 1
    df = (
        bf.overlap(subcomps.reset_index(), gtf1, how="inner", suffixes=("_", ""))
        .drop(columns=["chrom_", "start_", "end_"])
        .rename(columns=lambda c: c.rstrip("_"))
        .set_index("gene_id")
    )
    return de_table.merge(df, left_index=True, right_index=True)


def filter_by_gene_type(df: pd.DataFrame, *gene_types) -> pd.DataFrame:
    return df[df["gene_type"].str.lower().isin([s.lower() for s in gene_types])]


def compute_subcompartment_delta(df: pd.DataFrame, contrast: str, condition: str) -> npt.NDArray:
    return (df[condition].map(get_subcompartment_ranks()) - df[contrast].map(get_subcompartment_ranks())).to_numpy()


def compute_fisher_exact_test(df: pd.DataFrame, lfc_cutoff: float) -> Tuple[float, float]:
    df["delta"] = df["condition"].map(get_subcompartment_ranks()) - df["contrast"].map(get_subcompartment_ranks())

    m = np.zeros([2, 2], dtype=int)
    m[0][0] = len(df[(df["log2FoldChange"] >= lfc_cutoff) & (df["delta"] < 0)])
    m[0][1] = len(df[(df["log2FoldChange"] >= lfc_cutoff) & (df["delta"] > 0)])

    m[1][0] = len(df[(df["log2FoldChange"] < -lfc_cutoff) & (df["delta"] < 0)])
    m[1][1] = len(df[(df["log2FoldChange"] < -lfc_cutoff) & (df["delta"] > 0)])

    return ss.fisher_exact(m)


def plot_heatmap(df: pd.DataFrame, ax: plt.Axes, lfc_type: str, lfc_cutoff: float) -> npt.NDArray:
    if lfc_type == "down":
        df = df[df["log2FoldChange"] <= -lfc_cutoff]
    elif lfc_type == "up":
        df = df[df["log2FoldChange"] >= lfc_cutoff]
    elif lfc_type == "de":
        df = df[df["log2FoldChange"].abs() >= lfc_cutoff]
    else:
        assert lfc_type == "non-de"
        df = df[df["log2FoldChange"].abs() < lfc_cutoff]

    m = np.zeros([len(get_subcompartment_ranks()), len(get_subcompartment_ranks())], dtype=int)

    for (s1, s2), df1 in df.groupby(["contrast", "condition"]):
        if s1 == s2:
            continue
        i1 = get_subcompartment_ranks().get(s1)
        i2 = get_subcompartment_ranks().get(s2)
        m[i1, i2] = len(df1)

    img = ax.imshow(m)
    plt.colorbar(img, ax=ax)

    subcompartment_labels = list(get_subcompartment_ranks().keys())
    ax.set_xticks(range(len(subcompartment_labels)))
    ax.set_yticks(range(len(subcompartment_labels)))

    ax.set_xticklabels(subcompartment_labels)
    ax.set_yticklabels(subcompartment_labels)

    return m


def plot_diag_ratio(ax: plt.Axes, m: npt.NDArray):
    values = [1.0]
    for i in range(1, m.shape[0]):
        ratio = np.diag(m, i).sum() / np.diag(m, -i).sum()
        values.append(ratio)

    values = np.nan_to_num(values, nan=1.0, neginf=1.0, posinf=1.0)

    ax.plot(list(range(m.shape[0])), np.log2(values))
    ax.plot(list(range(m.shape[0])), [0] * len(values))
    ax.set(ylim=(-3.5, 3.5))
    ax.set_xticks(np.arange(0, len(values)))
    ax.set_xticklabels(np.arange(0, len(values)))


def plot_heatmaps(
    df: pd.DataFrame,
    contrast: str,
    condition: str,
    gene_types: List[str],
    lfc_cutoff: float,
    padj_cutoff: float,
):
    padj = "padj" if "padj" in df else "svalue"
    mask = df["gene_type"].str.lower().isin([s.lower() for s in gene_types])
    de = df[(df[padj] <= padj_cutoff) & mask].rename(columns={condition: "condition", contrast: "contrast"})
    nonde = df[(df[padj] > padj_cutoff) & mask].rename(columns={condition: "condition", contrast: "contrast"})

    types = ["down", "non-de", "up"]

    fig, axs = plt.subplots(2, len(types), figsize=(6.4 * len(types), 1.25 * 6.4), height_ratios=[1, 0.25])

    m1 = plot_heatmap(de, axs[0][0], "down", lfc_cutoff)
    m2 = plot_heatmap(nonde, axs[0][1], "non-de", lfc_cutoff)
    m3 = plot_heatmap(de, axs[0][2], "up", lfc_cutoff)

    plot_diag_ratio(axs[1][0], m1)
    plot_diag_ratio(axs[1][1], m2)
    plot_diag_ratio(axs[1][2], m3)

    for t, ax, m in zip(types, axs[0], (m1, m2, m3)):
        ax.set(title=f"{t} genes (sum={m.sum()}; tril={np.tril(m).sum()}; triu={np.triu(m).sum()})")

    df1 = df.rename(columns={condition: "condition", contrast: "contrast"})
    fisher_stat, fisher_pvalue = compute_fisher_exact_test(df1[df1[padj] <= padj_cutoff], lfc_cutoff)

    fig.suptitle(
        f"{contrast} vs {condition} (lfc={lfc_cutoff:.2f}, pval={padj_cutoff:.2f}, fisher_stat={fisher_stat:.2f}, fisher_pval={fisher_pvalue:.2g}): "
        + ";".join(gene_types)
    )
    fig.tight_layout()
    return fig


def build_pvalue_matrix_upreg(df: pd.DataFrame, lfc_cutoff: float) -> npt.NDArray:
    num_comps = len(get_subcompartment_ranks())
    m = np.full([num_comps, num_comps], np.nan, dtype=float)

    for s1, i1 in get_subcompartment_ranks().items():
        for s2, i2 in get_subcompartment_ranks().items():
            if i2 == i1:
                continue

            dff = df[(df["contrast"].isin({s1, s2})) & (df["condition"].isin({s1, s2}))]

            if i1 > i2:
                df1 = dff[(dff["log2FoldChange"] >= lfc_cutoff) & (dff["delta"] < 0)]
            else:
                df1 = dff[(dff["log2FoldChange"] >= lfc_cutoff) & (dff["delta"] > 0)]

            df2 = dff[(dff["log2FoldChange"] >= lfc_cutoff) & (dff["delta"] != 0)]

            if len(df2) == 0:
                continue

            res = ss.binomtest(len(df1), len(df2), alternative="greater")
            m[i1, i2] = res.pvalue

    return m


def build_pvalue_matrix_downreg(df: pd.DataFrame, lfc_cutoff: float) -> npt.NDArray:
    num_comps = len(get_subcompartment_ranks())
    m = np.full([num_comps, num_comps], np.nan, dtype=float)

    for s1, i1 in get_subcompartment_ranks().items():
        for s2, i2 in get_subcompartment_ranks().items():
            if i2 == i1:
                continue

            dff = df[(df["contrast"].isin({s1, s2})) & (df["condition"].isin({s1, s2}))]

            if i1 > i2:
                df1 = dff[(dff["log2FoldChange"] <= -lfc_cutoff) & (dff["delta"] < 0)]
            else:
                df1 = dff[(dff["log2FoldChange"] <= -lfc_cutoff) & (dff["delta"] > 0)]

            df2 = dff[(dff["log2FoldChange"] <= -lfc_cutoff) & (dff["delta"] != 0)]

            if len(df2) == 0:
                continue

            res = ss.binomtest(len(df1), len(df2), alternative="greater")
            m[i1, i2] = res.pvalue

    return m


def plot_pvalue_matrix(m: npt.NDArray, ax: plt.Axes):
    img = ax.imshow(m, vmin=0, vmax=0.05)

    for y in range(m.shape[0]):
        for x in range(m.shape[1]):
            ax.text(
                x,
                y,
                f"{m[y, x]:.2f}",
                ha="center",
                va="center",
            )

    subcompartment_labels = list(get_subcompartment_ranks().keys())
    ax.set_xticks(range(len(subcompartment_labels)))
    ax.set_yticks(range(len(subcompartment_labels)))

    ax.set_xticklabels(subcompartment_labels)
    ax.set_yticklabels(subcompartment_labels)

    plt.colorbar(img, ax=ax)


def plot_pvalue_heatmaps(
    df: pd.DataFrame,
    contrast: str,
    condition: str,
    gene_types: List[str],
    lfc_cutoff: float,
    padj_cutoff: float,
):
    padj = "padj" if "padj" in df else "svalue"
    mask = df["gene_type"].str.lower().isin([s.lower() for s in gene_types])
    de = df[(df[padj] <= padj_cutoff) & mask].rename(columns={condition: "condition", contrast: "contrast"})
    de["delta"] = de["condition"].map(get_subcompartment_ranks()) - de["contrast"].map(get_subcompartment_ranks())

    types = ["down", "up"]

    fig, axs = plt.subplots(1, len(types), figsize=(6.4 * len(types), 6.4))

    pv_downreg = build_pvalue_matrix_downreg(de, lfc_cutoff)
    pv_upreg = build_pvalue_matrix_upreg(de, lfc_cutoff)

    plot_pvalue_matrix(pv_downreg, axs[0])
    plot_pvalue_matrix(pv_upreg, axs[1])

    fig.suptitle(f"{contrast} vs {condition} (lfc={lfc_cutoff:.2f}, pval={padj_cutoff:.2f}): " + ";".join(gene_types))
    fig.tight_layout()
    return fig


def generate_line_collection(x, y, z, cmap="Blues_r"):
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/multicolored_line.html
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    colors = (np.array(z[1:]) + np.array(z[:-1])) / 2

    # norm = plt.Normalize(0.0, 0.1)
    norm = LogNorm(0.01, 1.0)
    lc = LineCollection(segments, cmap=cmap, norm=norm)
    lc.set_array(colors)

    return lc


def correlate(df: pd.DataFrame, key: str, x: npt.NDArray) -> Tuple[npt.NDArray, npt.NDArray, npt.NDArray, npt.NDArray]:
    pccs = []
    pvals = []
    nobs = []
    for cutoff in x:
        df1 = df[df[key].abs() >= cutoff]
        n = len(df1)
        if n > 1:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ss.ConstantInputWarning)
                pcc, pval = ss.pearsonr(df1["log2FoldChange"], df1["delta"])
            if np.isfinite(pcc):
                pccs.append(pcc)
                pvals.append(pval)
            else:
                pccs.append(0)
                pvals.append(1)
        else:
            pccs.append(0)
            pvals.append(1)
        nobs.append(n)

    mask = np.array(nobs) != 0
    return x[mask], np.array(pccs)[mask], np.array(pvals)[mask], np.array(nobs)[mask]


def plot_correlation_helper(x, pccs, pvals, nobs, fig, ax1):
    lc = generate_line_collection(x, pccs, pvals)
    line = ax1.add_collection(lc)
    fig.colorbar(line, ax=ax1)
    ax1.set_xlim(np.min(x), np.max(x))

    if np.min(pvals) > 0.1:
        ymin = -1.0
        ymax = 1.0
    else:
        ymin = min(0.0, pccs[pvals <= 0.1].min())
        ymax = pccs[pvals <= 0.1].max()

    ax1.set_ylim(ymin, ymax)

    ax2 = ax1.twinx()
    ax2.plot(x, nobs, color="orange")
    ax2.set(yscale="log", ylabel="# genes")


def plot_correlation(df: pd.DataFrame, contrast: str, condition: str, gene_types: List[str]) -> plt.Figure:
    df = filter_by_gene_type(df, *gene_types).copy()
    df["delta"] = compute_subcompartment_delta(df, contrast, condition)

    (fig, (ax11, ax21)) = plt.subplots(2, 1, figsize=(2 * 6.4, 6.4 * 1.5))

    x, pccs, pvals, nobs = correlate(df, "log2FoldChange", np.linspace(0, 15, 100))
    plot_correlation_helper(x, pccs, pvals, nobs, fig, ax11)
    ax11.set(ylabel="PCC", xlabel="log2FoldChange")

    x, pccs, pvals, nobs = correlate(df, "delta", np.arange(0, 8))
    plot_correlation_helper(x, pccs, pvals, nobs, fig, ax21)

    ax21.set(ylabel="PCC", xlabel="delta cutoff")
    fig.suptitle(f"{contrast} vs {condition}: " + ";".join(gene_types))

    return fig


def plot_boxplot(
    df: pd.DataFrame,
    contrast: str,
    condition: str,
    gene_types: List[str],
    lfc_cutoff: float,
    padj_cutoff: float,
) -> plt.Figure:
    df = filter_by_gene_type(df, *gene_types).copy()

    padj = "padj" if "padj" in df else "svalue"
    df = df[(df["log2FoldChange"].abs() >= lfc_cutoff) & (df[padj] < padj_cutoff)]

    df["delta"] = compute_subcompartment_delta(df, contrast, condition)
    max_delta = df["delta"].abs().max()

    colors = {}
    for i, n in zip(range(-max_delta, max_delta + 1), np.linspace(0.0, 1.0, max_delta * 2 + 1)):
        colors[i] = mpl.colormaps["bwr"](n)

    fig, ax = plt.subplots(1, 1)
    sns.boxplot(
        df[df["log2FoldChange"] < 0],
        x="delta",
        y="log2FoldChange",
        palette=colors,
        ax=ax,
    )
    sns.boxplot(
        df[df["log2FoldChange"] >= 0],
        x="delta",
        y="log2FoldChange",
        palette=colors,
        ax=ax,
    )

    fig.suptitle(f"{contrast} vs {condition} (padj={padj_cutoff:.2f}; lfc={lfc_cutoff:.2f}): " + ";".join(gene_types))

    return fig


def save_plot_to_file(fig: plt.Figure, outprefix: pathlib.Path, force: bool, close_after_save: bool = True) -> None:
    png = pathlib.Path(f"{outprefix}.png")
    svg = pathlib.Path(f"{outprefix}.svg")
    if not force:
        handle_path_collisions(png, svg)

    logging.info("Saving plot to %s.{png,svg}...", outprefix)
    fig.savefig(png, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight")
    if close_after_save:
        plt.close(fig)


def setup_logger(level=logging.INFO):
    fmt = "[%(asctime)s] %(levelname)s: %(message)s"
    logging.basicConfig(format=fmt)
    logging.getLogger().setLevel(level)


def main():
    args = vars(make_cli().parse_args())

    contrast = args["contrast"]
    condition = args["condition"]
    assert contrast != condition

    output_prefix = args["output-prefix"]
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    out_table = output_prefix.with_suffix(".tsv.gz")
    if not args["force"]:
        handle_path_collisions(out_table)

    df = (
        overlap_sucompartments_with_tss(
            import_subcomps(args["bedgraph"], contrast, condition),
            import_de_gene_table(args["de-table"], contrast, condition),
            import_gtf(args["gtf"]),
        )
        .rename(columns={"gene_name_y": "gene_name"})
        .drop(columns=["gene_name_x"])
    )

    logging.info(f"Writing table to {out_table}...")
    df.to_csv(out_table, sep="\t", index=False)

    if args["plot_type"] == "correlation":
        logging.info("Plotting correlation...")
        fig = plot_correlation(df, contrast, condition, args["gene_types"])
        save_plot_to_file(fig, output_prefix, args["force"])
    elif args["plot_type"] == "heatmap":
        logging.info("Plotting heatmap...")
        fig = plot_heatmaps(df, contrast, condition, args["gene_types"], args["lfc"], args["padj"])
        save_plot_to_file(fig, output_prefix, args["force"])
    elif args["plot_type"] == "pvalue-heatmap":
        logging.info("Plotting pvalue-heatmap...")
        fig = plot_pvalue_heatmaps(df, contrast, condition, args["gene_types"], args["lfc"], args["padj"])
        save_plot_to_file(fig, output_prefix, args["force"])
    elif args["plot_type"] == "boxplot":
        logging.info("Plotting boxplot...")
        fig = plot_boxplot(df, contrast, condition, args["gene_types"], args["lfc"], args["padj"])
        save_plot_to_file(fig, output_prefix, args["force"])
    else:
        assert False


if __name__ == "__main__":
    mpl.rcParams.update(
        {
            "axes.titlesize": 10,
            "axes.labelsize": 22,
            "legend.fontsize": 17,
            "xtick.labelsize": 18,
            "ytick.labelsize": 18,
        }
    )
    setup_logger()
    main()
