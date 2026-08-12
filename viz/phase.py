# /// script
# dependencies = ["numpy", "matplotlib"]
# ///
"""The wave function over spacetime, phase as hue and modulus as brightness.

The fronts of constant hue lean at the phase velocity, the bright band at the
group velocity, and for a free electron the first is half the second.
"""

import pathlib

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb

out = pathlib.Path(__file__).resolve().parent.parent / "out"
psi = np.load(out / "psi.npy")
position = np.load(out / "position.npy")
time = np.load(out / "time.npy")

order = position.argsort()
position = position[order]
psi = psi[:, order]

# The packet occupies a narrow band, and the fronts are only legible at the
# resolution the data itself has.
inside = (position > -32) & (position < 14)
position = position[inside]
psi = psi[:, inside]

hue = (np.angle(psi) + np.pi) / (2 * np.pi)
value = np.abs(psi) / np.abs(psi).max()
picture = hsv_to_rgb(np.stack([hue, np.ones_like(value), value], axis=-1))

figure, axes = plt.subplots(figsize=(7, 5), dpi=150, layout="constrained")
axes.imshow(
    picture,
    origin="lower",
    aspect="auto",
    extent=(position[0], position[-1], time[0], time[-1]),
    interpolation="nearest",
)
start = -30.0
for velocity, style, name in [(0.5788, "-", "group"), (0.2894, "--", "phase")]:
    axes.plot(
        start + velocity * time,
        time,
        style,
        color="white",
        linewidth=1.0,
        alpha=0.7,
        label=f"{name}, {velocity:.3f} nm/fs",
    )
axes.legend(loc="lower right", framealpha=0.3, labelcolor="white")
axes.set_xlim(position[0], position[-1])
axes.set_ylim(time[0], time[-1])
axes.set_xlabel("position [nm]")
axes.set_ylabel("time [fs]")

wheel = figure.colorbar(
    plt.cm.ScalarMappable(cmap="hsv", norm=plt.Normalize(-np.pi, np.pi)),
    ax=axes,
    ticks=[-np.pi, 0, np.pi],
)
wheel.ax.set_yticklabels([r"$-\pi$", "0", r"$\pi$"])
wheel.set_label(r"$\arg\psi$")

figure.savefig(out / "phase.png")
