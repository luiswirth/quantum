# /// script
# dependencies = ["numpy", "matplotlib", "pillow"]
# ///
"""The wave function as the curve it is, from the line into the complex plane.

A packet is a helix: its radius is the modulus, its pitch the wavelength, and
free evolution turns it while it drifts and stretches. The camera follows the
packet, so what is left to see is the turning.
"""

import pathlib

import matplotlib.pyplot as plt
import movie
import numpy as np

out = pathlib.Path(__file__).resolve().parent.parent / "out"
psi = np.load(out / "psi.npy")
position = np.load(out / "position.npy")
time = np.load(out / "time.npy")

order = position.argsort()
position = position[order]
psi = psi[:, order]


def refine(values, factor):
    """The same band limited function on a finer grid, by padding its
    spectrum, which invents nothing the samples did not already fix."""
    spectrum = np.fft.fft(values, axis=-1)
    half = values.shape[-1] // 2
    padded = np.zeros(values.shape[:-1] + (values.shape[-1] * factor,), complex)
    padded[..., :half] = spectrum[..., :half]
    padded[..., -half:] = spectrum[..., -half:]
    return np.fft.ifft(padded, axis=-1) * factor


refinement = 8
psi = refine(psi, refinement)
spacing = position[1] - position[0]
position = position[0] + spacing / refinement * np.arange(psi.shape[-1])

density = np.abs(psi) ** 2
center = (density * position).sum(1) / density.sum(1)
window = 8.0
height = np.abs(psi).max()

shape = (7.0, 4.2)
resolution = 150
size = movie.pixels(shape, resolution)
figure = plt.figure(figsize=shape, dpi=resolution)
axes = figure.add_axes((0.0, -0.04, 1.0, 1.0), projection="3d")
axes.view_init(elev=16, azim=-62)
axes.set_box_aspect((3.0, 1.0, 1.0), zoom=1.1)
axes.set_ylim(-1.15 * height, 1.15 * height)
axes.set_zlim(-1.15 * height, 1.15 * height)
# The complex plane is squeezed by the box aspect, so it carries three ticks
# and no more.
axes.set_yticks([-1.0, 0.0, 1.0])
axes.set_zticks([-1.0, 0.0, 1.0])
axes.tick_params(labelsize=8, pad=0)
axes.set_xlabel("position [nm]", labelpad=8)
axes.set_ylabel(r"$\Re\,\psi$", labelpad=2)
axes.set_zlabel(r"$\Im\,\psi$", labelpad=2)

(helix,) = axes.plot([], [], [], color="#1f4e79", linewidth=1.8)
(shadow,) = axes.plot([], [], [], color="#111111", linewidth=1.0, alpha=0.55)
(spine,) = axes.plot([], [], [], color="#999999", linewidth=0.8)
clock = figure.text(0.5, 0.93, "", ha="center", fontsize=12)

encoder = movie.encoder(out / "corkscrew.mp4", size)
for frame in range(len(time)):
    inside = np.abs(position - center[frame]) < window
    x = position[inside]
    value = psi[frame][inside]
    axes.set_xlim(center[frame] - window, center[frame] + window)
    helix.set_data_3d(x, value.real, value.imag)
    # The real part alone, laid on the far wall, is the picture a plot of
    # `Re psi` gives.
    shadow.set_data_3d(x, value.real, np.full_like(x, -1.15 * height))
    spine.set_data_3d(x, np.zeros_like(x), np.zeros_like(x))
    clock.set_text(f"t = {time[frame]:5.1f} fs")
    figure.canvas.draw()
    encoder.stdin.write(movie.capture(figure, size))

encoder.stdin.close()
encoder.wait()
figure.savefig(out / "corkscrew.png")
