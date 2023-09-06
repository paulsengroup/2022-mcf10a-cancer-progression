#!/usr/bin/env python3


# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import multiprocessing as mp
import pathlib
import subprocess as sp
import sys
import tempfile
import warnings
from typing import List

import bioframe as bf
import cooler
import cooltools
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import upsetplot


def make_cli():
    cli = argparse.ArgumentParser()

    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli.add_argument(
        "--loops",
        nargs=3,
        required=True,
        type=pathlib.Path,
        help="Path to three BEDPE with the list of loops to plot.",
    )

    cli.add_argument(
        "--mcools",
        nargs=3,
        required=True,
        type=str,
        help="Path to three cooler files (should be in the same order as the loops).",
    )

    cli.add_argument(
        "--lowest-resolution",
        required=True,
        type=positive_int,
        help="Lowest resolution used for loop calling.",
    )

    cli.add_argument(
        "--labels",
        type=str,
        nargs="+",
        help="List of labels to use instead of file names.",
    )

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=str,
        help="Output prefix.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )
    cli.add_argument(
        "--skip-header",
        action="store_true",
        default=False,
        help="Skip file header.",
    )
    cli.add_argument(
        "--nproc",
        type=int,
        choices=range(1, mp.cpu_count() + 1),
        default=mp.cpu_count(),
        help="Maximum number of parallel processes.",
    )

    cli.add_argument("--skip-pileup", action="store_true", default=False)

    return cli


def import_bedpe(path_to_bedpe: pathlib.Path, skip_header: bool) -> pd.DataFrame:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cols = bf.SCHEMAS["bedpe"][:6]
        return bf.read_table(path_to_bedpe, schema="bedpe", skiprows=1 if skip_header else 0)[cols]


def handle_path_collisions(*paths: pathlib.Path) -> None:
    collisions = [p for p in paths if p.exists()]

    if len(collisions) != 0:
        collisions = "\n - ".join((str(p) for p in collisions))
        raise RuntimeError(
            "Refusing to overwrite file(s):\n" f" - {collisions}\n" "Pass --force to overwrite existing file(s)."
        )


def save_plot_to_file(fig: plt.Figure, outprefix: str, force: bool, close_after_save: bool = True) -> None:
    outprefix = pathlib.Path(outprefix)
    png = outprefix.with_suffix(".png")
    svg = outprefix.with_suffix(".svg")
    if not force:
        handle_path_collisions(png, svg)

    fig.savefig(png, bbox_inches="tight", dpi=300)
    fig.savefig(svg, bbox_inches="tight")
    if close_after_save:
        plt.close(fig)


def pad_loops(df: pd.DataFrame, target_size: int) -> pd.DataFrame:
    df = df.copy()
    df["span"] = df["end1"] - df["start1"]
    dfs = []

    for span, grp in df.groupby("span"):
        if target_size < span:
            dfs.append(grp)
            continue

        padding = target_size - span
        grp = bf.expand(grp, padding // 2, cols=("chrom1", "start1", "end1"))
        grp = bf.expand(grp, padding // 2, cols=("chrom2", "start2", "end2"))

        dfs.append(grp.drop(columns="span"))

    return pd.concat(dfs)


def run_pair_to_pair(df1, df2, type="both", padding=0, frac=0.01):
    with tempfile.NamedTemporaryFile("wt") as bedpe1, tempfile.NamedTemporaryFile("wt") as bedpe2:
        df1.to_csv(bedpe1.name, sep="\t", header=False, index=False)
        df2.to_csv(bedpe2.name, sep="\t", header=False, index=False)

        cmd = [
            "pairToPair",
            "-a",
            bedpe1.name,
            "-b",
            bedpe2.name,
            "-is",
            "-type",
            type,
            "-f",
            str(frac),
        ]
        if padding != 0:
            cmd.extend(("-slop", str(padding)))
        with sp.Popen(cmd, stdout=sp.PIPE) as bedtools:
            df = pd.read_table(bedtools.stdout, usecols=list(range(6)), names=bf.SCHEMAS["bedpe"][:6])
            bedtools.communicate(timeout=60)

        if (code := bedtools.returncode) != 0:
            print(bedtools.stderr, file=sys.stderr)
            raise RuntimeError(f"bedtools terminated with code {code}")

        return df.drop_duplicates(keep="first")


def plot_upset_diagram(ab, bc, ac, a_only, b_only, c_only, shared, labels, output_prefix, force):
    subsets = pd.Series(
        {
            (True, False, False): len(a_only),
            (False, True, False): len(b_only),
            (False, False, True): len(c_only),
            (True, True, False): len(ab),
            (True, False, True): len(ac),
            (False, True, True): len(bc),
            (True, True, True): len(shared),
        }
    )
    subsets.index.names = labels
    fig, ax = plt.subplots()
    upset = upsetplot.UpSet(subsets)
    upset.plot(fig=fig)

    save_plot_to_file(fig, output_prefix, force)


def compute_pileup(df: pd.DataFrame, clr: cooler.Cooler, **kwargs) -> npt.NDArray:
    if "flank" not in kwargs:
        kwargs["flank"] = 100_000

    expected = compute_expected(clr, nproc=kwargs.get("nproc", 1))

    pileup = cooltools.api.snipping.pileup(
        clr,
        df,
        expected_df=expected,
        view_df=cooltools.lib.common.make_cooler_view(clr),
        **kwargs,
    )
    return np.nan_to_num(pileup).sum(axis=2)


def plot_pileups(
    output_prefix: str,
    clr: cooler.Cooler,
    a_only,
    b_only,
    c_only,
    ab,
    ac,
    bc,
    shared,
    labels: List[str],
    force: bool,
    nproc: int,
):
    def plot_pileup(pileup: npt.NDArray, title, output_prefix, force):
        fig, ax = plt.subplots(1, 1)
        ax.imshow(pileup, cmap="bwr")
        ax.set(title=title)
        save_plot_to_file(fig, output_prefix, force)

    plot_pileup(
        compute_pileup(shared, clr, nproc=nproc),
        "shared",
        output_prefix + "_shared",
        force,
    )
    plot_pileup(
        compute_pileup(ab, clr),
        f"{labels[0]}_{labels[1]}",
        output_prefix + f"_{labels[0]}_{labels[1]}",
        force,
    )
    plot_pileup(
        compute_pileup(bc, clr, nproc=nproc),
        f"{labels[1]}_{labels[2]}",
        output_prefix + f"_{labels[1]}_{labels[2]}",
        force,
    )
    plot_pileup(
        compute_pileup(ac, clr, nproc=nproc),
        f"{labels[0]}_{labels[2]}",
        output_prefix + f"_{labels[0]}_{labels[2]}",
        force,
    )

    plot_pileup(
        compute_pileup(a_only, clr, nproc=nproc),
        f"{labels[0]}_only",
        output_prefix + f"_{labels[0]}_only",
        force,
    )
    plot_pileup(
        compute_pileup(b_only, clr, nproc=nproc),
        f"{labels[1]}_only",
        output_prefix + f"_{labels[1]}_only",
        force,
    )
    plot_pileup(
        compute_pileup(c_only, clr, nproc=nproc),
        f"{labels[2]}_only",
        output_prefix + f"_{labels[2]}_only",
        force,
    )


@functools.cache
def compute_expected(clr, **kwargs) -> pd.DataFrame:
    return cooltools.api.expected.expected_cis(clr, **kwargs)


def main():
    args = vars(make_cli().parse_args())

    labels = args.get("labels")
    paths_to_loops = args["loops"]
    if len(paths_to_loops) != 3:
        raise RuntimeError(f"Expected path to 3 loop annotations, found {len(paths_to_loops)}")

    if labels is None:
        labels = [p.name for p in paths_to_loops]
    elif len(labels) != len(paths_to_loops):
        raise RuntimeError(
            f"Mismatch in the number of files and labels: expected {len(paths_to_loops)} labels, found {len(labels)}"
        )

    loops = {
        k: pad_loops(import_bedpe(path, args["skip_header"]), args["lowest_resolution"])
        for k, path in zip(labels, paths_to_loops)
    }

    base_res = min(
        cooler.Cooler(args["mcools"][0] + "::" + suffix).binsize
        for suffix in cooler.fileops.list_coolers(args["mcools"][0])
    )

    clrs = {k: cooler.Cooler(p + f"::/resolutions/{base_res}") for k, p in zip(labels, args["mcools"])}

    clra = clrs[labels[0]]
    clrb = clrs[labels[1]]
    clrc = clrs[labels[2]]

    a = loops[labels[0]]
    b = loops[labels[1]]
    c = loops[labels[2]]

    shared = run_pair_to_pair(run_pair_to_pair(a, b), run_pair_to_pair(a, c))

    ab = run_pair_to_pair(a, b)
    ab = run_pair_to_pair(ab, shared, "notboth")

    ac = run_pair_to_pair(a, c)
    ac = run_pair_to_pair(ac, shared, "notboth")

    bc = run_pair_to_pair(b, c)
    bc = run_pair_to_pair(bc, shared, "notboth")

    a_only = run_pair_to_pair(a, b, "notboth")
    a_only = run_pair_to_pair(a_only, c, "notboth")

    b_only = run_pair_to_pair(b, a, "notboth")
    b_only = run_pair_to_pair(b_only, c, "notboth")

    c_only = run_pair_to_pair(c, a, "notboth")
    c_only = run_pair_to_pair(c_only, b, "notboth")

    plot_upset_diagram(
        ab,
        bc,
        ac,
        a_only,
        b_only,
        c_only,
        shared,
        labels,
        args["output_prefix"] + "_upset",
        args["force"],
    )

    if args["skip_pileup"]:
        return

    for label, clr in zip(labels, (clra, clrb, clrc)):
        plot_pileups(
            args["output_prefix"] + f"_{label}",
            clr,
            a_only,
            b_only,
            c_only,
            ab,
            ac,
            bc,
            shared,
            labels,
            args["force"],
            args["nproc"],
        )


if __name__ == "__main__":
    main()
