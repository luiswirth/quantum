# /// script
# dependencies = ["numpy", "matplotlib"]
# ///
"""The probability density over spacetime, time upward and space across."""

import pathlib

import matplotlib.pyplot as plt
import numpy as np

out = pathlib.Path(__file__).resolve().parent.parent / "out"
psi = np.load(out / "psi.npy")
position = np.load(out / "position.npy")
time = np.load(out / "time.npy")

order = position.argsort()
density = np.abs(psi[:, order]) ** 2

figure, axes = plt.subplots(figsize=(6, 5), layout="constrained")
axes.imshow(
    density,
    origin="lower",
    aspect="auto",
    cmap="magma",
    extent=(position[order][0], position[order][-1], time[0], time[-1]),
)
axes.set_xlabel("position [nm]")
axes.set_ylabel("time [fs]")
figure.savefig(out / "spacetime.png", dpi=160)
