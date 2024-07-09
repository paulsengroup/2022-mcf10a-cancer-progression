#!/usr/bin/env python3

# Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

import argparse
import logging
import os.path
import pathlib
from typing import Any, Dict, List, Tuple

import cv2
import h5py
import numpy as np
import numpy.typing as npt
import pandas as pd
import PIL
import scipy
import skimage
from PIL.Image import Image


def make_cli():
    def positive_int(arg):
        if (n := int(arg)) > 0:
            return n

        raise ValueError("Not a positive integer")

    cli = argparse.ArgumentParser()

    cli.add_argument(
        "figures",
        nargs="+",
        type=pathlib.Path,
        help="Path to one or more figures to be processed.",
    )

    cli.add_argument(
        "-o",
        "--output",
        type=pathlib.Path,
        required=True,
        help="Path where to save the processed images in HDF5 format.",
    )

    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Overwrite existing files (if any).",
    )

    cli.add_argument(
        "--crop-nuclei",
        action="store_true",
        default=False,
        help="Crop nuclei images removing empty spaces in the image.",
    )

    cli.add_argument(
        "--contour-padding",
        type=positive_int,
        default=50,
        help="Padding used to expand contours of individual nuclei.",
    )

    cli.add_argument(
        "--sharpen-radius",
        type=float,
        default=5.0,
        help="Sharpen radius.\n"
        "See https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.unsharp_mask for more details.",
    )

    cli.add_argument(
        "--sharpen-amount",
        type=float,
        default=5.0,
        help="Sharpen amount.\n"
        "See https://scikit-image.org/docs/stable/api/skimage.filters.html#skimage.filters.unsharp_mask for more details.",
    )

    cli.add_argument(
        "--min-object-size",
        type=int,
        default=5000,
        help="Cutoff on the object size used to discard contours of small objects.",
    )

    return cli


def open_image(path: str, scale: float = 1.0) -> Image:
    img = PIL.Image.open(path)

    if scale == 1:
        return img

    size = int((float(img.size[1]) * scale))
    return img.resize((size, size), PIL.Image.Resampling.LANCZOS)


def trim_zeros(arr):
    """
    https://stackoverflow.com/questions/55917328/numpy-trim-zeros-in-2d-or-3d
    Returns a trimmed view of an n-D array excluding any outer
    regions which contain only zeros.
    """
    slices = tuple(slice(idx.min(), idx.max() + 1) for idx in np.nonzero(arr))
    return arr[slices]


def segment_nuclei(
    img: Image,
    sharpen_radius: float,
    sharpen_amount: float,
    min_obj_size: int,
    contour_padding: int,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # Make grayscale
    _, _, blue = img.split()
    blue = np.array(blue)

    # Sharpen
    blue = skimage.filters.unsharp_mask(blue, radius=sharpen_radius, amount=sharpen_amount)

    # Find objects in foreground
    threshold = skimage.filters.threshold_otsu(blue)
    mask = blue > threshold
    mask = skimage.morphology.remove_small_objects(mask, min_obj_size)
    labels = skimage.measure.label(mask)

    props = skimage.measure.regionprops(labels, blue)

    # Compute objects contours
    contours = []
    for index in range(0, labels.max()):
        label_i = props[index].label
        contour = skimage.measure.find_contours(labels == label_i, 0.5)[0]

        y, x = contour.T
        contours.append(pd.DataFrame({"x": x, "y": y, "label": [label_i] * len(x)}))

    # Expand labels and re-compute object contours
    labels = skimage.segmentation.expand_labels(labels, distance=contour_padding)
    props = skimage.measure.regionprops(labels, blue)

    # Compute objects contours
    contours_expanded = []
    for index in range(0, labels.max()):
        label_i = props[index].label
        contour = skimage.measure.find_contours(labels == label_i, 0.5)[0]

        y, x = contour.T
        contours_expanded.append(pd.DataFrame({"x": x, "y": y, "label": [label_i] * len(x)}))

    return pd.concat(contours), pd.concat(contours_expanded)


def crop_image_based_on_contours(img: Image, contours: pd.DataFrame, crop_nuclei: bool) -> List[npt.NDArray]:
    imgs = []

    for label, df in contours.groupby("label"):
        # Create an empty image to store the masked array
        mask = np.zeros_like(img.split()[0], dtype="bool")

        # Create a contour image by using the contour coordinates rounded to their nearest integer value
        x = df["x"].round().astype("int")
        y = df["y"].round().astype("int")

        # Ensure contour is closed
        if x.iloc[0] != x.iloc[-1] or y.iloc[0] != y.iloc[-1]:
            continue

        mask[y, x] = True

        # Fill in the hole created by the contour boundary
        mask = scipy.ndimage.binary_fill_holes(mask)

        np_img = np.array(img, dtype=np.uint8)
        np_img[~mask] = 0

        if crop_nuclei:
            np_img = trim_zeros(np_img)

        imgs.append(np_img)

    return imgs


def score_nuclei(img: npt.NDArray, contours: pd.DataFrame) -> npt.NDArray:
    scores = []

    for _, df in contours.groupby("label"):
        if len(df) < 6:
            scores.append(np.nan)
            continue

        contour_coo = df[["x", "y"]].rename(columns={"x": "y", "y": "x"}).to_numpy().astype(int)
        min_ellipse = cv2.fitEllipse(contour_coo)

        ellipse_2d = np.zeros(img.shape, dtype=np.uint8)
        contour_2d = np.zeros(img.shape, dtype=np.uint8)

        cv2.drawContours(contour_2d, [contour_coo], 0, 1, thickness=-1)
        cv2.ellipse(ellipse_2d, min_ellipse, 1, thickness=-1)

        ellipse_2d = ellipse_2d.sum(axis=2)
        contour_2d = contour_2d.sum(axis=2)

        score = (ellipse_2d & contour_2d).sum() / min(ellipse_2d.sum(), contour_2d.sum())
        scores.append(score)

    return np.array(scores, dtype=float)


def process_image(path: str, args: Dict[str, Any]) -> Tuple[npt.NDArray, List[npt.NDArray], pd.DataFrame]:
    img = open_image(path)
    contours, contours_expanded = segment_nuclei(
        img,
        sharpen_radius=args["sharpen_radius"],
        sharpen_amount=args["sharpen_amount"],
        min_obj_size=args["min_object_size"],
        contour_padding=args["contour_padding"],
    )

    imgs = crop_image_based_on_contours(img, contours_expanded, args["crop_nuclei"])

    contours_ = []
    for (_, df1), (_, df2) in zip(contours.groupby("label"), contours_expanded.groupby("label")):
        # Ensure contour is closed
        if df2["x"].iloc[0] != df2["x"].iloc[-1] or df2["y"].iloc[0] != df2["y"].iloc[-1]:
            continue

        df1["label"] = len(contours_)
        df1["x"] -= df2["x"].min()
        df1["y"] -= df2["y"].min()

        contours_.append(df1)

    if len(contours_) == 0:
        contours_ = [pd.DataFrame([], columns=["label", "x", "y"], dtype=int)]

    return np.array(img), [np.array(img_) for img_ in imgs], pd.concat(contours_)


def write_array_to_hdf5(m: npt.NDArray, h5: h5py.File, path: str):
    h5.create_dataset(path, data=m, compression="gzip", compression_opts=9)


def main():
    args = vars(make_cli().parse_args())

    output = args["output"]
    if not args["force"] and os.path.exists(output):
        raise RuntimeError(f'Refusing to overwrite output file "{output}"')

    with h5py.File(output, "w") as h5:
        for p in args["figures"]:
            logging.info(f"Processing {os.path.basename(p)}")
            try:
                input_img, imgs, contours = process_image(p, args)
            except ValueError as e:
                if "not enough values to unpack" in str(e):
                    logging.warning(f"Image {os.path.basename(p)} is not an RGB image. SKIPPING!")
                    continue
                raise
            if len(contours) == 0:
                continue

            scores = score_nuclei(input_img, contours)

            prefix = str(p).replace("\\", "_").replace("/", "_").lstrip("_")
            write_array_to_hdf5(input_img, h5, f"{prefix}/imgs/000_input")

            for (
                i,
                img,
            ) in enumerate(imgs):
                write_array_to_hdf5(img, h5, f"{prefix}/imgs/{i+1:03d}")

            write_array_to_hdf5(scores, h5, f"{prefix}/scores")

            write_array_to_hdf5(contours["x"], h5, f"{prefix}/contours/x")
            write_array_to_hdf5(contours["y"], h5, f"{prefix}/contours/y")
            write_array_to_hdf5(contours["label"], h5, f"{prefix}/contours/label")


def setup_logger(level=logging.INFO):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s: %(message)s")
    logging.getLogger().setLevel(level)


if __name__ == "__main__":
    setup_logger()
    main()
