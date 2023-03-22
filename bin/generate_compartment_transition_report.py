#!/usr/bin/env python3

# Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import functools
import itertools
import pathlib
import subprocess as sp
import sys
import tempfile
from collections import Counter
from typing import Dict, Tuple, Union

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from numpy.lib.stride_tricks import sliding_window_view


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    def positive_float(s) -> float:
        if (x := float(s)) > 0:
            return x

        raise RuntimeError("Not a positive float")

    cli.add_argument(
        "bedgraph",
        type=pathlib.Path,
        help="Path to the (sub)compartment bedGraph generated by dcHiC.",
    )

    cli.add_argument(
        "--plot-type",
        nargs="+",
        default=["alluvial", "heatmap"],
        type=str,
        choices={"alluvial", "heatmap"},
        required=True,
        help="Type of plot to generate.",
    )

    cli.add_argument(
        "--path-to-plotting-script",
        type=pathlib.Path,
        default="make_ab_comp_alluvial.r",
    )

    cli.add_argument(
        "-o",
        "--output-prefix",
        type=pathlib.Path,
        required=True,
        help="Path to output prefix (including parent folder(s) but without extension).",
    )

    cli.add_argument(
        "--aggregate-subcompartments",
        action="store_true",
        default=False,
        help="Aggregate subcompartments into A/B compartments.",
    )

    cli.add_argument(
        "--base-color",
        type=str,
        default="grey",
        help="Default color used in alluvial plots.",
    )

    cli.add_argument(
        "--highlight-color",
        type=str,
        help="Color used for highlighting in alluvial plots.",
    )

    cli.add_argument(
        "--width",
        type=positive_float,
        default=6.4,
        help="Plot width.",
    )

    cli.add_argument(
        "--height",
        type=positive_float,
        default=6.4,
        help="Plot height.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwrite existing files.",
    )

    return cli


@functools.cache
def get_compartment_ranks() -> dict:
    compartment_labels = tuple(["B", "B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3", "A"])
    return {k: v for v, k in enumerate(compartment_labels)}


@functools.cache
def get_compartment_color_mappings() -> dict:
    compartment_labels = tuple(["B", "B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3", "A"])
    colors = tuple(
        [
            "#3b4cc0ff",
            "#3b4cc0ff",
            "#6788eeff",
            "#9abbffff",
            "#c9d7f0ff",
            "#edd1c2ff",
            "#f7a889ff",
            "#e26952ff",
            "#b40426ff",
            "#b40426ff",
        ]
    )
    return {k: v for k, v in zip(compartment_labels, colors)}


@functools.cache
def get_comp_to_label_mappings() -> dict:
    ranks = get_compartment_ranks()
    highest_rank = max(ranks.values())
    padding_size = len(str(highest_rank))

    return {comp: f"__{i:0{padding_size}d}__" for comp, i in ranks.items()}


@functools.cache
def get_label_to_comp_mappings() -> dict:
    return {v: k for k, v in get_comp_to_label_mappings().items()}


def import_data(path_to_df: pathlib.Path) -> Tuple[int, pd.DataFrame]:
    if str(path_to_df) == "-":
        path_to_df = sys.stdin

    df = pd.read_table(path_to_df).rename(columns={"chr": "chrom"})

    bin_size = (df["end"] - df["start"]).max()
    df.loc[df["state.mode"] == "None", "state.mode"] = pd.NA
    return bin_size, df.set_index(["chrom", "start", "end"])


def group_and_sort_subcompartments(df: pd.DataFrame, aggregate_subcompartments: bool) -> pd.DataFrame:
    if aggregate_subcompartments:
        df = df.apply(lambda x: x.str[0])

    df = df.groupby(df.columns.tolist(), as_index=False).size()

    cols = [col for col in df.columns if col != "size"]
    df.sort_values(
        by=cols,
        key=lambda x: x.apply(lambda y: get_compartment_ranks().get(y)),
        ignore_index=True,
        inplace=True,
        kind="stable",
    )

    return df


def get_compartment_labels_from_df(df: pd.DataFrame) -> list:
    return pd.concat([df[col] for col in df.columns if col != "size"]).unique().tolist()


def rename_compartments(df: pd.DataFrame, mappings: Union[dict, None] = None) -> pd.DataFrame:
    """
    We have to rename compartments such that similar compartments will be next to each other in the alluvial plot
    The plotting library in R allows to customize stream ordering, but this feature seems very brittle.
    """
    df = df.copy()
    if mappings is None:
        mappings = get_comp_to_label_mappings()

    for col in df.columns:
        if col.endswith(".state"):
            df[col] = df[col].apply(lambda x: mappings[x])

    return df


def plot_alluvial_with_r(
    path_to_input: pathlib.Path,
    output_name: pathlib.Path,
    path_to_rscript: pathlib.Path,
    base_color: str,
    highlight_color: str,
    width: float,
    height: float,
    compartment_to_highlight: Union[str, None],
) -> None:
    cmd = [
        str(path_to_rscript),
        f"--compartments={path_to_input}",
        f"--base_color={base_color}",
        f"--highlight_color={highlight_color}",
        f"--width={width}",
        f"--height={height}",
        "--outprefix",
        str(output_name.with_suffix("")),
    ]
    if compartment_to_highlight is not None:
        cmd.append(f"--highlight_label={compartment_to_highlight}")

    if base_color == "white":
        cmd.append("--only_show_highlighted")

    try:
        sp.run(cmd, capture_output=True, text=True).check_returncode()
    except sp.CalledProcessError as e:
        cmd = " ".join((x for x in cmd))
        print(
            f'Process "{cmd}" terminated with exit code {e.returncode}.',
            file=sys.stderr,
        )
        if len(e.stdout) != 0:
            print(f"STDOUT:\n{e.stdout}", file=sys.stderr)
        if len(e.stderr) != 0:
            print(f"STDERR:\n{e.stderr}", file=sys.stderr)

        raise e


def replace_compartment_labels_svg(path_to_svg: pathlib.Path):
    data = path_to_svg.read_text()
    assert len(data) != 0

    for label, comp in get_label_to_comp_mappings().items():
        data = data.replace(label, comp)

    with path_to_svg.open("w") as f:
        f.write(data)


def make_alluvial_plot(
    path_to_tsv: pathlib.Path,
    compartment_to_highlight: str,
    outname: pathlib.Path,
    path_to_plotting_script: pathlib.Path,
    overwrite_existing: bool,
    base_color: Union[str, None],
    highlight_color: Union[str, None],
    width: float,
    height: float,
) -> None:
    if outname.exists() and not overwrite_existing:
        raise RuntimeError(f"Refusing to overwrite file {outname}. Pass --force to overwrite existing file(s).")

    comp_label = get_comp_to_label_mappings().get(compartment_to_highlight)

    if base_color is None:
        base_color = "white"
    if highlight_color is None:
        if comp_label is not None:
            highlight_color = get_compartment_color_mappings()[get_label_to_comp_mappings().get(comp_label)]
        else:
            highlight_color = "orange"

    plot_alluvial_with_r(
        path_to_input=path_to_tsv,
        output_name=outname,
        path_to_rscript=path_to_plotting_script,
        base_color=base_color,
        highlight_color=highlight_color,
        width=width,
        height=height,
        compartment_to_highlight=comp_label,
    )

    replace_compartment_labels_svg(outname)


def compute_subcomp_transision_matrix_cis(
    subcompartment_labels: pd.Series, subcomp_ranks: Dict[str, int]
) -> npt.NDArray:
    counter = Counter()
    for _, labels in subcompartment_labels.groupby(level="chrom"):
        pairs = sliding_window_view(labels, 2)
        counter += Counter(zip(pairs[:, 0], pairs[:, 1]))

    counts = np.zeros([len(subcomp_ranks), len(subcomp_ranks)], dtype=int)
    for (sc1, sc2), count in counter.items():
        i1 = subcomp_ranks[sc1]
        i2 = subcomp_ranks[sc2]
        counts[i1, i2] = count

    return counts


def compute_subcomp_transition_matrix_trans(
    states1: pd.Series,
    states2: pd.Series,
    subcomp_ranks: Dict[str, int],
) -> np.ndarray:
    assert len(states1) == len(states2)

    counts = np.zeros([len(subcomp_ranks), len(subcomp_ranks)], dtype=int)
    for (sc1, sc2), count in Counter(zip(states1, states2)).items():
        i1 = subcomp_ranks[sc1]
        i2 = subcomp_ranks[sc2]
        counts[i1, i2] = count

    return counts


def make_heatmap_helper(grid: np.ndarray, ax: plt.Axes, cmap: str, label_scale_factor: float = 1.0):
    # Negative values will be transparent
    cmap = mpl.cm.get_cmap(cmap)
    cmap.set_under("k", alpha=0)

    ax.imshow(grid, cmap=cmap, vmin=0)

    # Show Mbps involved in transitions
    max_val = np.nanmax(grid)
    for (j, i), label in np.ndenumerate(grid):
        if label < 0:
            continue
        if label / max_val >= 0.40:
            color = "white"
        else:
            color = "black"

        # I can't figure out a better way to make sure a float is formatted using at most 4 digits (plus the dot)
        label = str(label / label_scale_factor)
        if len(label) > 5:
            label = label[:5]
        ax.text(i, j, label, ha="center", va="center", color=color)


def make_heatmap_trans(
    states1: pd.Series,
    states2: pd.Series,
    bin_size: int,
    ax: plt.Axes,
    cmap1: str = "Reds",
    cmap2: str = "Greys",
) -> None:
    if len(set(states1.unique().tolist() + states2.unique().tolist())) > 2:
        ranks = {}
        for comp in get_compartment_ranks().keys():
            if comp != "A" and comp != "B":
                ranks[comp] = len(ranks)
    else:
        ranks = {"B": 0, "A": 1}

    label1 = states1.name
    label2 = states2.name
    assert label1 != label2

    grid = compute_subcomp_transition_matrix_trans(states1, states2, ranks)
    grid *= bin_size

    rsum = np.sum(grid, axis=0)
    csum = np.sum(grid, axis=1)

    # Plot grid
    filler = np.full_like(grid[0], -1.0)
    grid1 = np.c_[grid, filler, filler]  # Add two transparent columns to the right
    make_heatmap_helper(grid1, ax, cmap1, 1.0e6)

    # Plot transition frequencies (normalized by rowsum + colsum)
    transition_vect = (2 * np.diag(grid)) / (rsum + csum)
    # Make all but the sencond to last column transparent
    filler1 = np.full_like(grid, -1.0)
    filler2 = np.full_like(transition_vect, -1.0)
    grid2 = np.c_[filler1, transition_vect, filler2]
    make_heatmap_helper(grid2, ax, cmap2)

    # Plot transition frequencies (normalized by matrix sum)
    transition_vect = np.diag(grid) / grid.sum()
    # Make all but the last column transparent
    filler1 = np.full_like(grid, -1.0)
    filler2 = np.full_like(transition_vect, -1.0)
    grid3 = np.c_[filler1, filler2, transition_vect]
    make_heatmap_helper(grid3, ax, cmap2)

    ax.set(title=f"Compartment transitions (Mbp) - {label1} vs {label2}", xlabel=label2, ylabel=label1)

    ax.set_xticks(range(len(ranks)), labels=list(ranks.keys()))
    ax.set_yticks(range(len(ranks)), labels=list(ranks.keys()))


def make_heatmap_cis(
    states: pd.Series,
    bin_size: int,
    ax: plt.Axes,
    cmap1: str = "Reds",
    cmap2: str = "Greys",
) -> None:
    if states.nunique() > 2:
        ranks = {}
        for comp in get_compartment_ranks().keys():
            if comp != "A" and comp != "B":
                ranks[comp] = len(ranks)
    else:
        ranks = {"B": 0, "A": 1}

    label = states.name
    grid = compute_subcomp_transision_matrix_cis(states, ranks)
    grid *= bin_size

    rsum = np.sum(grid, axis=0)
    csum = np.sum(grid, axis=1)

    # Plot diagonal only
    grid1 = np.full_like(grid, -1.0)
    np.fill_diagonal(grid1, grid.diagonal())
    filler = np.full_like(grid[0], -1.0)
    grid1 = np.c_[grid1, filler]  # Add one transparent column to the right
    make_heatmap_helper(grid1, ax, cmap2, 1.0e6)

    # Plot grid except diagonal
    filler = np.full_like(grid[0], -1.0)
    grid2 = grid.copy()
    np.fill_diagonal(grid2, -1.0)
    grid2 = np.c_[grid2, filler]  # Add one transparent columns to the right
    make_heatmap_helper(grid2, ax, cmap1, 1.0e6)

    # Plot transition frequencies (normalized by row/col)
    transition_vect = (2 * np.diag(grid)) / (rsum + csum)
    filler = np.full_like(grid, -1.0)
    grid3 = np.c_[filler, transition_vect]
    make_heatmap_helper(grid3, ax, cmap2)

    ax.set(title=f"Sub-compartment transitions (Mbp) - {label}", xlabel=label, ylabel=label)

    ax.set_xticks(range(len(ranks)), labels=list(ranks.keys()))
    ax.set_yticks(range(len(ranks)), labels=list(ranks.keys()))


def compute_transition_coefficients(df: pd.DataFrame) -> pd.DataFrame:
    labels = get_compartment_labels_from_df(df)

    data = {"label": [], "same": [], "switch": []}
    for label in labels:
        df1 = df[(df == label).any(axis="columns")]
        mask = (df1.drop(columns=["size"]) == label).all(axis="columns")
        data["label"].append(label)
        data["same"].append(df1.loc[mask, "size"].sum())
        data["switch"].append(df1.loc[~mask, "size"].sum())

    df1 = pd.DataFrame(data)
    df1["ratio"] = df1["same"] / (df1["same"] + df1["switch"])

    return df1


def handle_file_name_collision(outname: pathlib.Path, force: bool) -> None:
    if outname.exists() and not force:
        raise RuntimeError(f"Refusing to overwrite file {outname}. Pass --force to overwrite existing file(s).")


def add_path_suffix(path: pathlib.Path, suffix: str, extension: str) -> pathlib.Path:
    return pathlib.Path(path.parent) / f"{path.name}{suffix}{extension}"


def make_alluvial_plots_subcmd(args: Dict) -> None:
    _, df = import_data(args["bedgraph"])
    df = group_and_sort_subcompartments(df.filter(regex=r"\.state$"), args["aggregate_subcompartments"])

    output_prefix = args["output_prefix"]

    compartment_labels = get_compartment_labels_from_df(df)
    with tempfile.NamedTemporaryFile() as tmpdata:
        rename_compartments(df).to_csv(tmpdata.name, sep="\t", index=False)

        for comp in compartment_labels:
            outname = output_prefix.parent / f"{output_prefix.name}_{comp}.svg"
            make_alluvial_plot(
                path_to_tsv=pathlib.Path(tmpdata.name),
                compartment_to_highlight=comp,
                outname=outname,
                path_to_plotting_script=args["path_to_plotting_script"],
                overwrite_existing=args["force"],
                base_color=args.get("base_color"),
                highlight_color=args.get("highlight_color"),
                width=args["width"],
                height=args["height"],
            )


def make_heatmap_plots_subcmd(args: Dict) -> None:
    bin_size, df = import_data(args["bedgraph"])

    if args["aggregate_subcompartments"]:
        conditions = df.filter(regex=".state$").columns
        df[conditions] = df[conditions].apply(lambda x: x.str[0])

    output_prefix = args["output_prefix"]

    conditions = df.filter(regex=".state$").columns.tolist()
    for cond1, cond2 in itertools.combinations_with_replacement(conditions, 2):
        fig, ax = plt.subplots(1, 1, figsize=(args["height"], args["width"]))
        if cond1 == cond2:
            make_heatmap_cis(df[cond1], bin_size, ax)
            cond1 = cond1.removesuffix(".state")
            outstem = output_prefix.parent / f"{output_prefix.name}_{cond1}_subcomp_transition_heatmap"
        else:
            make_heatmap_trans(df[cond1], df[cond2], bin_size, ax)
            cond1 = cond1.removesuffix(".state")
            cond2 = cond2.removesuffix(".state")
            outstem = output_prefix.parent / f"{output_prefix.name}_{cond1}_vs_{cond2}_subcomp_transition_heatmap"

        if (outstem.with_suffix(".png").exists() or outstem.with_suffix(".svg").exists()) and not args["force"]:
            raise RuntimeError(
                f"Refusing to overwrite file: {outstem}.{{png,svg}}\n" "Pass --force to overwrite existing file(s)."
            )

        plt.tight_layout()
        fig.savefig(outstem.with_suffix(".svg"))
        fig.savefig(outstem.with_suffix(".png"), dpi=600)


def main():
    args = vars(make_cli().parse_args())
    args["output_prefix"].parent.mkdir(parents=True, exist_ok=True)

    if "alluvial" in args["plot_type"]:
        make_alluvial_plots_subcmd(args)
    if "heatmap" in args["plot_type"]:
        make_heatmap_plots_subcmd(args)


if __name__ == "__main__":
    main()
