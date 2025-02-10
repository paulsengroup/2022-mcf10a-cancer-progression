#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import itertools
import os
import pathlib
import sys
import warnings
from typing import Any, Dict, List, Optional, Union

import cv2
import h5py
import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages


def make_cli():
    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "URI",
        type=str,
        help="HDF5 URI to a figure generated by one of the scripts in this repository.",
    )

    cli.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        required=True,
        help="Path where to save the extracted image.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    cli.add_argument(
        "--plot-blobs",
        action="store_true",
        default=False,
        help="Overlay blobs on top of the original image.",
    )

    cli.add_argument(
        "--channel",
        choices={"RGB", "R", "G", "B"},
        default="RGB",
        help="Image channels to be plotted.",
    )

    cli.add_argument(
        "--overlap-cfx-cutoff",
        default=0.0,
        type=float,
        help="Overlap coefficient cutoff used to filter out nuclei whose shape diverges from that of an ellipsis.",
    )

    cli.add_argument(
        "--dpi",
        type=positive_int,
        default=300,
        help="Figure DPI.",
    )

    cli.add_argument(
        "--colorblind-friendly",
        action="store_true",
        default=False,
        help="Re-paint image using a more color-blind friendly color map.",
    )

    cli.add_argument(
        "--red-mapping",
        nargs=3,
        type=int,
        help="RGB code used for mapping the R channel of the given image.",
    )
    cli.add_argument(
        "--green-mapping",
        nargs=3,
        type=int,
        help="RGB code used for mapping the G channel of the given image.",
    )
    cli.add_argument(
        "--blue-mapping",
        nargs=3,
        type=int,
        help="RGB code used for mapping the B channel of the given image.",
    )
    cli.add_argument(
        "--black-cutoff",
        type=int,
        help="Cutoff used to detect black pixels.",
    )
    cli.add_argument(
        "--white-cutoff",
        type=int,
        help="Cutoff used to detect white pixels.",
    )

    return cli


def import_blobs(grp: Union[h5py.Group, None], color: str) -> Union[pd.DataFrame, None]:
    if grp is None:
        return None

    return pd.DataFrame({"x": grp[f"{color}/x"][:], "y": grp[f"{color}/y"][:]})


def plot_many(
    images_h5: List[h5py.Dataset],
    blobs_h5: List[Union[h5py.Group, None]],
    channel: str,
    plot_blobs: bool,
    red_mapping: Optional[List[int]],
    green_mapping: Optional[List[int]],
    blue_mapping: Optional[List[int]],
    black_cutoff: Optional[int],
    white_cutoff: Optional[int],
) -> plt.Figure:
    cols = 4
    rows = max(1, int(np.ceil(len(images_h5) / cols)))
    fig1, axs = plt.subplots(rows, cols, figsize=(6.4 * cols, 6.4 * rows))

    try:
        _ = next(itertools.chain(*axs))
    except:
        axs = [axs]

    for ax in itertools.chain(*axs):
        ax.axis("off")

    if len(images_h5) != 1:
        blobs_h5 = [None] + blobs_h5

    for img_dset, blobs_grp, ax in zip(images_h5, blobs_h5, itertools.chain(*axs)):
        if img_dset is None:
            print("Skipping image...", file=sys.stderr)
            continue
        img = recolor_image(img_dset[:], red_mapping, green_mapping, blue_mapping, black_cutoff, white_cutoff)
        red_blobs = import_blobs(blobs_grp, color="red")
        green_blobs = import_blobs(blobs_grp, color="green") if plot_blobs else None
        if channel == "R":
            print("Plotting RED channel...", file=sys.stderr)
            ax.imshow(img[:, :, 0], cmap="gray")
        elif channel == "G":
            print("Plotting GREEN channel...", file=sys.stderr)
            ax.imshow(img[:, :, 1], cmap="gray")
        elif channel == "B":
            print("Plotting BLUE channel...", file=sys.stderr)
            ax.imshow(img[:, :, 2], cmap="gray")
        else:
            print("Plotting RGB image...", file=sys.stderr)
            ax.imshow(img)

        if red_blobs is not None:
            if red_mapping is None:
                color = "red"
            else:
                color = np.array(red_mapping, dtype=float)
                color /= color.sum()
            print("Plotting red blobs...", file=sys.stderr)
            ax.scatter(
                red_blobs["x"],
                red_blobs["y"],
                s=20,
                color=color,
                edgecolors="white",
                linewidth=1,
            )

        if green_blobs is not None:
            if green_mapping is None:
                color = "green"
            else:
                color = np.array(green_mapping, dtype=float)
                color /= color.sum()
            print("Plotting green blobs...", file=sys.stderr)
            ax.scatter(
                green_blobs["x"],
                green_blobs["y"],
                s=20,
                color=color,
                edgecolors="white",
                linewidth=1,
            )

    fig1.tight_layout()
    return fig1


def process_one_image(uri: str, args: Dict[str, Any]) -> plt.Figure:
    file_path, _, hdf5_path = uri.partition("::")
    hdf5_img_path = f"{hdf5_path}/imgs/"
    hdf5_blobs_path = f"{hdf5_path}/blobs/"

    with h5py.File(file_path) as h5:
        obj = h5.get(hdf5_img_path)
        is_dataset = isinstance(obj, h5py.Dataset)

        if obj is None:
            raise RuntimeError(f"Object does not exist: {hdf5_img_path}")

        if is_dataset:
            img_objects = [obj]
            blobs_objects = [h5.get(hdf5_blobs_path) if args["plot_blobs"] else None]
        else:
            assert isinstance(obj, h5py.Group)
            img_objects = [h5[f"{hdf5_img_path}/{grp}"] for grp in obj]
            if args["plot_blobs"]:
                blobs_objects = [h5.get(f"{hdf5_blobs_path}/{dset}") for dset in h5.get(hdf5_blobs_path)]
            else:
                blobs_objects = [None] * len(img_objects)

            scores = h5[f"{hdf5_path}/scores"][:]

            for i, score in enumerate(scores):
                if score < args["overlap_cfx_cutoff"]:
                    img_objects[i + 1] = None
                    blobs_objects[i] = None

        fig = plot_many(
            img_objects,
            blobs_objects,
            args["channel"],
            args["plot_blobs"],
            red_mapping=args["red_mapping"],
            green_mapping=args["green_mapping"],
            blue_mapping=args["blue_mapping"],
            black_cutoff=args["black_cutoff"],
            white_cutoff=args["white_cutoff"],
        )

        return fig


def recolor_image(
    img: npt.NDArray,
    red_mapping: Optional[List[int]],
    green_mapping: Optional[List[int]],
    blue_mapping: Optional[List[int]],
    black_cutoff: Optional[int],
    white_cutoff: Optional[int],
):
    if red_mapping is None:
        red_mapping = [255, 0, 0]
    if green_mapping is None:
        green_mapping = [0, 255, 0]
    if blue_mapping is None:
        blue_mapping = [0, 0, 255]
    if black_cutoff is None:
        black_cutoff = 0
    if white_cutoff is None:
        white_cutoff = 255

    red_mapping = np.array(red_mapping, dtype=float)
    green_mapping = np.array(green_mapping, dtype=float)
    blue_mapping = np.array(blue_mapping, dtype=float)

    red_mapping /= red_mapping.sum()
    green_mapping /= green_mapping.sum()
    blue_mapping /= blue_mapping.sum()

    black_cutoff /= 255
    white_cutoff /= 255

    red_channel, green_channel, blue_channel = cv2.split(img / 255)

    red_channel_3d = np.dstack([red_channel] * 3)
    green_channel_3d = np.dstack([green_channel] * 3)
    blue_channel_3d = np.dstack([blue_channel] * 3)

    red_channel = red_channel_3d[:, :, 0].T
    green_channel = green_channel_3d[:, :, 0].T
    blue_channel = blue_channel_3d[:, :, 0].T

    white_mask = (red_channel > white_cutoff) & (green_channel > white_cutoff) & (blue_channel > white_cutoff)
    black_mask = (red_channel < black_cutoff) & (green_channel < black_cutoff) & (blue_channel < black_cutoff)

    new_red_channel, new_green_channel, new_blue_channel = (
        (red_mapping * red_channel_3d) + (green_mapping * green_channel_3d) + (blue_mapping * blue_channel_3d)
    ).T

    new_red_channel[white_mask | black_mask] = red_channel[white_mask | black_mask]
    new_green_channel[white_mask | black_mask] = green_channel[white_mask | black_mask]
    new_blue_channel[white_mask | black_mask] = blue_channel[white_mask | black_mask]

    img = cv2.merge((new_blue_channel.T, new_green_channel.T, new_red_channel.T))

    return cv2.cvtColor(img.astype(np.float32), cv2.COLOR_BGR2RGB)


def main():
    args = vars(make_cli().parse_args())

    output = args["output"]
    if not args["force"] and os.path.exists(output):
        raise RuntimeError(f'Refusing to overwrite output file "{output}"')

    file_path, _, hdf5_path = args["URI"].partition("::")

    if args["colorblind_friendly"]:
        if args["red_mapping"] is not None:
            warnings.warn("--red-mapping is ignored when --colorblind-friendly is specified")
        args["red_mapping"] = [7, 82, 238]
        if args["green_mapping"] is not None:
            warnings.warn("--green-mapping is ignored when --colorblind-friendly is specified")
        args["green_mapping"] = [247, 247, 0]
        if args["blue_mapping"] is not None:
            warnings.warn("--blue-mapping is ignored when --colorblind-friendly is specified")
        args["blue_mapping"] = [245, 249, 249]
        if args["black_cutoff"] is not None:
            warnings.warn("--black-cutoff is ignored when --colorblind-friendly is specified")
        if args["white_cutoff"] is not None:
            warnings.warn("--white-cutoff is ignored when --colorblind-friendly is specified")

    if hdf5_path != "":
        fig = process_one_image(args["URI"], args)
        fig.tight_layout()
        fig.savefig(args["output"], dpi=args["dpi"])
        return

    if str(args["output"]).split(".")[-1] != "pdf":
        raise RuntimeError("Output format should be PDF when plotting multiple images at once.")

    with PdfPages(args["output"]) as pdf:
        with h5py.File(file_path) as h5:
            for grp in h5:
                fig = process_one_image(f"{file_path}::/{grp}", args)
                fig.suptitle(grp)
                fig.tight_layout()
                pdf.savefig(fig, dpi=args["dpi"])
                plt.close(fig)


if __name__ == "__main__":
    main()
