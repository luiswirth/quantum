"""Matplotlib frames piped into an h264 file."""

import subprocess

import numpy as np
from PIL import Image


def pixels(shape, resolution):
    """The even pixel count the encoder needs."""
    return tuple(2 * round(length * resolution / 2) for length in shape)


def encoder(path, size, frames_per_second=25):
    command = (
        "ffmpeg -y -loglevel error"
        f" -f rawvideo -pix_fmt rgb24 -s {size[0]}x{size[1]}"
        f" -framerate {frames_per_second} -i -"
        " -c:v libx264 -preset slow -crf 15 -pix_fmt yuv420p -movflags +faststart"
    ).split() + [str(path)]
    return subprocess.Popen(command, stdin=subprocess.PIPE)


def capture(figure, size):
    """The canvas is rendered at the display's pixel ratio, so the picture is
    brought to the size asked for before it is handed over."""
    picture = Image.fromarray(np.asarray(figure.canvas.buffer_rgba())).convert("RGB")
    return picture.resize(size, Image.LANCZOS).tobytes()
