#!/usr/bin/env python3

# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


import argparse
import pathlib

try:
    import cv2
    import dash_daq as daq
    import numpy as np
    import plotly.graph_objects as go
    from dash import Dash, Input, Output, dcc, html
except ImportError as e:
    raise RuntimeError(
        "Please make sure you have installed the following dependencies with pip: opencv-python, numpy, dash_daq, dash, plotly, kaleido."
    ) from e

app = Dash(__name__)


def make_cli() -> argparse.ArgumentParser:
    cli = argparse.ArgumentParser()

    cli.add_argument(
        "image",
        type=pathlib.Path,
        help="Path to the image to be processed.",
    )
    cli.add_argument(
        "--output",
        type=pathlib.Path,
        help="Path where to store the processed plot.",
    )
    cli.add_argument(
        "--force",
        action="store_true",
        default=False,
        help="Force overwriting of output file.",
    )
    cli.add_argument(
        "--red-mapping",
        nargs=3,
        type=int,
        default=[7, 82, 238],
        help="RGB code used for mapping the R channel of the given image.",
    )
    cli.add_argument(
        "--green-mapping",
        nargs=3,
        type=int,
        default=[247, 247, 0],
        help="RGB code used for mapping the G channel of the given image.",
    )
    cli.add_argument(
        "--blue-mapping",
        nargs=3,
        type=int,
        default=[245, 249, 249],
        help="RGB code used for mapping the B channel of the given image.",
    )
    cli.add_argument(
        "--black-cutoff",
        type=int,
        default=20,
        help="Cutoff used to detect black pixels.",
    )
    cli.add_argument(
        "--white-cutoff",
        type=int,
        default=220,
        help="Cutoff used to detect white pixels.",
    )

    grp = cli.add_mutually_exclusive_group()
    grp.add_argument(
        "--output-red-channel",
        action="store_true",
        default=False,
        help="Output the red channel only as a grayscale image.",
    )
    grp.add_argument(
        "--output-green-channel",
        action="store_true",
        default=False,
        help="Output the green channel only as a grayscale image.",
    )
    grp.add_argument(
        "--output-blue-channel",
        action="store_true",
        default=False,
        help="Output the blue channel only as a grayscale image.",
    )

    return cli


def mapping_slow(red_mapping, green_mapping, blue_mapping, black_cutoff, white_cutoff):
    assert black_cutoff <= white_cutoff

    new_red_channel = np.full_like(red_channel, 1.0)
    new_green_channel = np.full_like(green_channel, 1.0)
    new_blue_channel = np.full_like(blue_channel, 1.0)

    for i in range(red_channel.shape[0]):
        for j in range(red_channel.shape[1]):
            r = red_channel[i, j]
            g = green_channel[i, j]
            b = blue_channel[i, j]

            r, g, b = red_mapping * r + green_mapping * g + blue_mapping * b
            new_red_channel[i, j] = r
            new_green_channel[i, j] = g
            new_blue_channel[i, j] = b

    img = cv2.merge((new_blue_channel, new_green_channel, new_red_channel))
    return cv2.cvtColor(img, cv2.COLOR_BGR2RGB)


def mapping_fast(red_mapping, green_mapping, blue_mapping, black_cutoff, white_cutoff):
    assert black_cutoff <= white_cutoff

    red_channel = red_channel_3d[:, :, 0].T
    green_channel = green_channel_3d[:, :, 0].T
    blue_channel = blue_channel_3d[:, :, 0].T

    if white_cutoff > 1:
        white_cutoff /= 255
    if black_cutoff > 1:
        black_cutoff /= 255

    white_mask = (red_channel > white_cutoff) & (green_channel > white_cutoff) & (blue_channel > white_cutoff)
    black_mask = (red_channel < black_cutoff) & (green_channel < black_cutoff) & (blue_channel < black_cutoff)

    new_red_channel, new_green_channel, new_blue_channel = (
        (red_mapping * red_channel_3d) + (green_mapping * green_channel_3d) + (blue_mapping * blue_channel_3d)
    ).T

    new_red_channel[white_mask | black_mask] = red_channel[white_mask | black_mask]
    new_green_channel[white_mask | black_mask] = green_channel[white_mask | black_mask]
    new_blue_channel[white_mask | black_mask] = blue_channel[white_mask | black_mask]

    img = cv2.merge((new_blue_channel.T, new_green_channel.T, new_red_channel.T)) * 255

    return cv2.cvtColor(img.astype(np.float32), cv2.COLOR_BGR2RGB)


def rgb_dict_to_numpy(rgb):
    x = np.array(
        [
            rgb["r"],
            rgb["g"],
            rgb["b"],
        ],
        dtype=float,
    )

    x /= x.sum()
    return x


def rgb_tuple_to_dict(rgb):
    assert len(rgb) == 3

    return dict(r=rgb[0], g=rgb[1], b=rgb[2], a=0)


def make_colorblind_friendly(
    red_mapping,
    green_mapping,
    blue_mapping,
    black_cutoff,
    white_cutoff,
):
    red_mapping = rgb_dict_to_numpy(red_mapping)
    green_mapping = rgb_dict_to_numpy(green_mapping)
    blue_mapping = rgb_dict_to_numpy(blue_mapping)

    return mapping_fast(
        red_mapping,
        green_mapping,
        blue_mapping,
        black_cutoff,
        white_cutoff,
    )


def setup_server(red_mapping, green_mapping, blue_mapping, black_cutoff, white_cutoff):
    app.layout = html.Div(
        [
            dcc.Graph(id="image-graph"),
            html.Div(
                [
                    daq.ColorPicker(
                        id="color-picker-red",
                        label="Red replacement",
                        value=dict(rgb=rgb_tuple_to_dict(red_mapping)),
                        size=164,
                    ),
                    daq.ColorPicker(
                        id="color-picker-green",
                        label="Green replacement",
                        value=dict(rgb=rgb_tuple_to_dict(green_mapping)),
                        size=164,
                    ),
                    daq.ColorPicker(
                        id="color-picker-blue",
                        label="Blue replacement",
                        value=dict(rgb=rgb_tuple_to_dict(blue_mapping)),
                        size=164,
                    ),
                ],
                style={"display": "flex", "justify-content": "space-around"},
            ),
            html.Div(
                [
                    html.H2("Interval Selection Slider"),
                    html.P(
                        "Use the slider below to select a range of values between 0 and 255.\n"
                        "Pixels whose RGB values are all outside the given range are not processed.\n"
                        "This is useful to keep blacks black, and whites white after mapping colors."
                    ),
                    dcc.RangeSlider(
                        0,
                        255,
                        1,
                        value=[black_cutoff, white_cutoff],
                        marks={i: str(i) for i in range(0, 256, 25)},
                        id="cutoff-slider",
                    ),
                ]
            ),
        ]
    )
    app.run_server(debug=True)


@app.callback(
    Output("image-graph", "figure"),
    [
        Input("color-picker-red", "value"),
        Input("color-picker-green", "value"),
        Input("color-picker-blue", "value"),
        Input("cutoff-slider", "value"),
    ],
)
def update_output(red_mapping, green_mapping, blue_mapping, cutoff_slider):
    img_mapped = make_colorblind_friendly(
        red_mapping=red_mapping["rgb"],
        green_mapping=green_mapping["rgb"],
        blue_mapping=blue_mapping["rgb"],
        black_cutoff=cutoff_slider[0],
        white_cutoff=cutoff_slider[1],
    )
    fig = go.Figure()
    fig.add_trace(go.Image(z=img_mapped))

    fig.update_layout(xaxis=dict(visible=False), yaxis=dict(visible=False), margin=dict(l=0, r=0, t=0, b=0))

    return fig


def load_image_raster(path):
    return cv2.imread(path).astype(np.float32) / 255.0


def load_image_svg(path):
    try:
        import cairosvg

        png_data = cairosvg.svg2png(url=str(path))
        img_array = np.frombuffer(png_data, np.uint8)
        img = cv2.imdecode(img_array, cv2.IMREAD_COLOR)
        return cv2.cvtColor(img, cv2.COLOR_BGRA2BGR) / 255.0
    except ImportError as e:
        raise RuntimeError("Processing images in svg format requires cairosvg to be installed") from e


def load_image_pdf(path):
    try:
        from pdf2image import convert_from_path

        img = np.array(convert_from_path(path)[0])
        return cv2.cvtColor(img, cv2.COLOR_RGBA2BGR) / 255.0
    except ImportError as e:
        raise RuntimeError("Processing images in pdf format requires pdf2image") from e


def load_image(path):
    assert path.exists()
    if path.suffix == ".pdf":
        return load_image_pdf(path)
    if path.suffix == ".svg":
        return load_image_svg(path)
    return load_image_raster(path)


if __name__ == "__main__":
    args = vars(make_cli().parse_args())

    if args["black_cutoff"] > args["white_cutoff"]:
        raise RuntimeError("--black-cutoff cannot be greater than --white-cutoff")

    img = load_image(args["image"])
    blue_channel, green_channel, red_channel = cv2.split(img)

    red_channel_3d = np.dstack([red_channel] * 3)
    green_channel_3d = np.dstack([green_channel] * 3)
    blue_channel_3d = np.dstack([blue_channel] * 3)

    dest = args["output"]

    if dest is None:
        setup_server(
            args["red_mapping"],
            args["green_mapping"],
            args["blue_mapping"],
            black_cutoff=args["black_cutoff"],
            white_cutoff=args["white_cutoff"],
        )
    else:
        if dest.exists() and not args["force"]:
            raise RuntimeError(f'Refusing to overwrite existing file "{dest}". Pass --force to overwrite.')

        img = make_colorblind_friendly(
            rgb_tuple_to_dict(args["red_mapping"]),
            rgb_tuple_to_dict(args["green_mapping"]),
            rgb_tuple_to_dict(args["blue_mapping"]),
            black_cutoff=args["black_cutoff"],
            white_cutoff=args["white_cutoff"],
        )

        if args["output_red_channel"]:
            img = img[:, :, 0]
        elif args["output_green_channel"]:
            img = img[:, :, 1]
        elif args["output_blue_channel"]:
            img = img[:, :, 2]

        fig = go.Figure()

        if len(img.shape) == 2:
            fig.add_trace(go.Heatmap(z=img, colorscale="gray"))
            fig.update_yaxes(autorange="reversed", scaleanchor="x", constrain="domain")
            fig.update_xaxes(constrain="domain")
        else:
            fig.add_trace(go.Image(z=img))
        fig.update_layout(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            margin=dict(l=0, r=0, t=0, b=0),
        )

        fig.write_image(dest)
