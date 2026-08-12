# /// script
# dependencies = ["numpy", "matplotlib", "pillow"]
# ///
"""The wave function itself, the camera following the packet.

The envelope travels at the group velocity and the ripples inside it at the
phase velocity, so for a free electron the ripples fall backwards through their
own envelope at half the speed.
"""

import pathlib

import matplotlib.pyplot as plt
import movie
import numpy as np
from matplotlib.colors import hsv_to_rgb
from matplotlib.patches import Polygon

out = pathlib.Path(__file__).resolve().parent.parent / "out"
psi = np.load(out / "psi.npy")
position = np.load(out / "position.npy")
time = np.load(out / "time.npy")

order = position.argsort()
position = position[order]
psi = psi[:, order]

density = np.abs(psi) ** 2
center = (density * position).sum(1) / density.sum(1)
window = 12.0
height = np.abs(psi).max()

# A fixed axes rectangle, since an automatic layout would resize it whenever
# the tick labels change width.
shape = (6.5, 3.6)
resolution = 150
size = movie.pixels(shape, resolution)
figure = plt.figure(figsize=shape, dpi=resolution)
axes = figure.add_axes((0.10, 0.14, 0.87, 0.78))
# The phase paints the body of the packet, and the modulus bounds it.
hue = axes.imshow(
    np.zeros((1, len(position), 3)),
    origin="lower",
    aspect="auto",
    extent=(position[0], position[-1], -1.15 * height, 1.15 * height),
    interpolation="nearest",
    zorder=0,
)
envelope = Polygon(np.zeros((2, 2)), closed=True, visible=False)
axes.add_patch(envelope)
hue.set_clip_path(envelope)
(real,) = axes.plot(
    [], [], color="black", linewidth=2.0, zorder=3, label=r"$\Re\,\psi$"
)
(upper,) = axes.plot(
    [], [], color="#111111", linewidth=1.6, zorder=3, label=r"$\pm|\psi|$"
)
(lower,) = axes.plot([], [], color="#111111", linewidth=1.6, zorder=3)
axes.set_ylim(-1.15 * height, 1.15 * height)
axes.set_xlabel("position [nm]")
axes.legend(loc="upper right")
clock = axes.set_title("")

encoder = movie.encoder(out / "packet.mp4", size)

for frame in range(len(time)):
    axes.set_xlim(center[frame] - window, center[frame] + window)
    real.set_data(position, psi[frame].real)
    upper.set_data(position, np.abs(psi[frame]))
    lower.set_data(position, -np.abs(psi[frame]))
    angle = (np.angle(psi[frame]) + np.pi) / (2 * np.pi)
    hue.set_data(
        hsv_to_rgb(np.stack([angle, np.ones_like(angle), np.ones_like(angle)], -1))[
            None
        ]
    )
    modulus = np.abs(psi[frame])
    envelope.set_xy(
        np.concatenate(
            [
                np.column_stack([position, modulus]),
                np.column_stack([position[::-1], -modulus[::-1]]),
            ]
        )
    )
    hue.set_clip_path(envelope)
    clock.set_text(f"t = {time[frame]:5.1f} fs")
    figure.canvas.draw()
    encoder.stdin.write(movie.capture(figure, size))

encoder.stdin.close()
encoder.wait()
figure.savefig(out / "packet.png")
