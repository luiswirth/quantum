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
import numpy as np
from PIL import Image

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
figure = plt.figure(figsize=(6.5, 3.6), dpi=150)
axes = figure.add_axes((0.10, 0.14, 0.87, 0.78))
(real,) = axes.plot([], [], color="#3b6ea5", linewidth=1.2, label=r"$\Re\,\psi$")
(upper,) = axes.plot([], [], color="#111111", linewidth=1.6, label=r"$\pm|\psi|$")
(lower,) = axes.plot([], [], color="#111111", linewidth=1.6)
axes.set_ylim(-1.15 * height, 1.15 * height)
axes.set_xlabel("position [nm]")
axes.legend(loc="upper right")
clock = axes.set_title("")


def draw(frame):
    axes.set_xlim(center[frame] - window, center[frame] + window)
    real.set_data(position, psi[frame].real)
    upper.set_data(position, np.abs(psi[frame]))
    lower.set_data(position, -np.abs(psi[frame]))
    clock.set_text(f"t = {time[frame]:5.1f} fs")
    figure.canvas.draw()
    return Image.frombuffer(
        "RGBA", figure.canvas.get_width_height(), figure.canvas.buffer_rgba()
    ).convert("RGB")


frames = [draw(index) for index in range(len(time))]
frames[0].save(
    out / "packet.webp",
    save_all=True,
    append_images=frames[1:],
    duration=40,
    loop=0,
    lossless=True,
    method=6,
)
frames[len(frames) // 2].save(out / "packet.png")
