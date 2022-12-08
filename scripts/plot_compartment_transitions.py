#!/usr/bin/env python3

import argparse
import functools
import pathlib
import subprocess as sp
import sys
import tempfile
from typing import Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "bedgraph",
        type=pathlib.Path,
        help="Path to the (sub)compartment bedGraph generated by dcHiC.",
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
        "--no-plotting",
        default=False,
        action="store_true",
        help="Skip plot generation. Only output a TSV with the data after aggregation.",
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
        help="Default color.",
    )

    cli.add_argument(
        "--highlight-color",
        type=str,
        default="orange",
        help="Color used for highlighting.",
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
    compartment_labels = tuple(
        ["B", "B3", "B2", "B1", "B0", "A0", "A1", "A2", "A3", "A"]
    )
    return {k: v for v, k in enumerate(compartment_labels)}


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

    df = pd.read_table(path_to_df)
    bin_size = (df["end"] - df["start"]).max()

    return bin_size, df.filter(regex=".state$")


def group_and_sort_subcompartments(
    df: pd.DataFrame, aggregate_subcompartments: bool
) -> pd.DataFrame:
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


def rename_compartments(
    df: pd.DataFrame, mappings: Union[dict, None] = None
) -> pd.DataFrame:
    """
    We have to rename compartments such that similar compartments will be next to each other in the alluvial plot
    The plotting library in R allows to customize stream ordering, but this feature seems very brittle.
    """
    df = df.copy()
    if mappings is None:
        mappings = get_comp_to_label_mappings()

    cols = [col for col in df.columns if col != "size"]
    for col in cols:
        df[col] = df[col].apply(lambda x: mappings[x])

    return df


def plot_alluvial_with_r(
    path_to_input: pathlib.Path,
    output_name: pathlib.Path,
    path_to_rscript: pathlib.Path,
    base_color: str,
    highlight_color: str,
    compartment_to_highlight: Union[str, None] = None,
) -> None:
    cmd = [
        str(path_to_rscript),
        "--compartments",
        str(path_to_input),
        "--base_color",
        base_color,
        "--highlight_color",
        highlight_color,
        "--outprefix",
        str(output_name.with_suffix("")),
    ]
    if compartment_to_highlight is not None:
        cmd.append(f"--highlight_label={compartment_to_highlight}")

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
    base_color: str,
    highlight_color: str,
) -> None:
    if outname.exists() and not overwrite_existing:
        raise RuntimeError(
            f"Refusing to overwrite file {outname}. Pass --force to overwrite existing file(s)."
        )

    comp_label = get_comp_to_label_mappings().get(compartment_to_highlight)
    plot_alluvial_with_r(
        path_to_input=path_to_tsv,
        output_name=outname,
        path_to_rscript=path_to_plotting_script,
        base_color=base_color,
        highlight_color=highlight_color,
        compartment_to_highlight=comp_label,
    )

    replace_compartment_labels_svg(outname)


def make_heatmap(
    df: pd.DataFrame,
    condition1: str,
    condition2: str,
    bin_size: int,
    ax: plt.Axes,
    cmap: str = "Reds",
) -> None:
    assert condition1 in df.columns
    assert condition2 in df.columns

    if len(df[condition1].unique()) > 2:
        ranks = {}
        for comp in get_compartment_ranks().keys():
            if comp != "A" and comp != "B":
                ranks[comp] = len(ranks)
    else:
        ranks = {"B": 0, "A": 1}

    grid = np.zeros([len(ranks), len(ranks)], dtype=int)

    for ((comp1, comp2), size) in (
        df.groupby([condition1, condition2])["size"].sum().items()
    ):
        i1, i2 = ranks[comp1], ranks[comp2]
        grid[i1, i2] = size * bin_size

    ax.imshow(grid, cmap=cmap)

    min_val = grid.min()
    delta = grid.max() - min_val

    for (j, i), label in np.ndenumerate(grid):
        if label - min_val >= (delta * 0.85):
            color = "white"
        else:
            color = "black"
        ax.text(i, j, f"{(label / 1.0e6):.1f}", ha="center", va="center", color=color)

    ax.set(title=f"Compartment transitions (Mbp) - {condition1} vs {condition2}")

    ax.set_xticks(range(len(ranks)), labels=list(ranks.keys()))
    ax.set_yticks(range(len(ranks)), labels=list(ranks.keys()))


def main():
    args = vars(make_cli().parse_args())

    bin_size, df = import_data(args["bedgraph"])
    df = group_and_sort_subcompartments(df, args["aggregate_subcompartments"])

    output_prefix = pathlib.Path(args["output_prefix"])
    output_prefix.parent.mkdir(exist_ok=True)

    df.to_csv(output_prefix.with_suffix(".tsv"), sep="\t", index=False)
    if args["no_plotting"]:
        return

    compartment_labels = get_compartment_labels_from_df(df)
    with tempfile.NamedTemporaryFile() as tmpdata:
        df1 = rename_compartments(df)
        df1.to_csv(tmpdata.name, sep="\t", index=False)

        for comp in compartment_labels:
            outname = output_prefix.parent / f"{output_prefix.name}_{comp}.svg"
            make_alluvial_plot(
                path_to_tsv=pathlib.Path(tmpdata.name),
                compartment_to_highlight=comp,
                outname=outname,
                path_to_plotting_script=args["path_to_plotting_script"],
                overwrite_existing=args["force"],
                base_color=args["base_color"],
                highlight_color=args["highlight_color"],
            )

    conditions = [col for col in df.columns if col != "size"]
    cond_pairs = []
    for i, cond1 in enumerate(conditions):
        for j, cond2 in enumerate(conditions):
            if i < j:
                cond_pairs.append(tuple([cond1, cond2]))

    num_plots = len(cond_pairs)
    fig, axs = plt.subplots(num_plots, 1, figsize=(6.4, 6.4 * num_plots))

    for (ax, (cond1, cond2)) in zip(axs, cond_pairs):
        make_heatmap(df, cond1, cond2, bin_size, ax)

    outname = (
        output_prefix.parent
        / f"{output_prefix.name}_compartment_transition_heatmaps.svg"
    )
    fig.savefig(outname)


if __name__ == "__main__":
    main()
