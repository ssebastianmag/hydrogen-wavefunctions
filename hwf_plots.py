"""
File: hwf_plots.py
Description: Plotting utilities for hydrogenic wavefunctions.
"""

from datetime import datetime
from typing import Optional, Literal
from pydantic import BaseModel

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.patheffects as pe
import seaborn as sns
import numpy as np

import hydrogen_wavefunction as hwf


class WaveFunction(BaseModel):
    """ Wavefunction schema.
        [Composite of 'hydrogen_wavefunctions.compute_psi_xz_slice'].
    """
    n: int
    l: int
    m: int
    Z: int = 1
    use_reduced_mass: bool = True
    M: Optional[float] = None
    extent_a_mu: float = 20.0
    grid_points: int = 600
    phi_value: float = 0.0
    phi_mode: Literal["plane", "constant"] = "plane"


def plot_hydrogen_wavefunction_xz(
        wf: WaveFunction,
        colormap: str = "rocket",
        use_dark_theme: bool = False,
        k: Optional[float] = None
):
    """ Plot hydrogen wavefunction restricted to the y=0 (x–z) plane.

        Parameters:
            wf (WaveFunction): Wavefunction parameters.
            colormap (str): Seaborn colormap name.
            use_dark_theme (bool): Plot theme rendering mode.
            k (float): Framing scale factor for extent calculation.
    """
    try:
        _ = sns.color_palette(colormap)
    except Exception:
        raise ValueError(f"(!) {colormap!r} is not a recognized Seaborn colormap.")

    if k is not None:  # Compute extent from framing scale factor (if provided)
        wf.extent_a_mu = float(k * (3 * wf.n * wf.n - wf.l * (wf.l + 1)) / (2 * wf.Z))

    # Retrieve X-Z grid, psi, reduced Bohr radius and P = |psi|^2
    Xg, Zg, psi, a_mu = hwf.compute_psi_xz_slice(
        n=wf.n, l=wf.l, m=wf.m, Z=wf.Z,
        use_reduced_mass=wf.use_reduced_mass,
        M=wf.M,
        extent_a_mu=wf.extent_a_mu,
        grid_points=wf.grid_points,
        phi_value=wf.phi_value,
        phi_mode=wf.phi_mode,
    )

    P = hwf.compute_probability_density(psi)

    # Global styles
    plt.rcParams["font.family"] = "STIXGeneral"
    plt.rcParams["mathtext.fontset"] = "stix"

    plt.rcParams["xtick.major.width"] = 4
    plt.rcParams["ytick.major.width"] = 4

    plt.rcParams["xtick.major.size"] = 15
    plt.rcParams["ytick.major.size"] = 15

    plt.rcParams["xtick.labelsize"] = 30
    plt.rcParams["ytick.labelsize"] = 30
    plt.rcParams["axes.linewidth"] = 4

    fig, ax = plt.subplots(figsize=(17, 16.5))
    plt.subplots_adjust(top=0.82, right=0.87, left=-0.10)

    # Theme styles
    if use_dark_theme:

        # Background color -> Darkest color in the colormap
        pal_100 = sns.color_palette(colormap, n_colors=100)
        background_color = sorted(pal_100, key=lambda c: 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2])[0]

        title_color = text_color = "#dfdfdf"
        tick_color = "#c4c4c4"

        fig.patch.set_facecolor(background_color)
        for spine in ax.spines.values():
            spine.set_color(tick_color)

        ax.tick_params(axis="x", colors=tick_color)
        ax.tick_params(axis="y", colors=tick_color)

    else:
        title_color = text_color = tick_color = "#000000"
        ax.tick_params(axis="x", colors=tick_color)
        ax.tick_params(axis="y", colors=tick_color)

    cmap = sns.color_palette(colormap, as_cmap=True)

    # Render plot within bounds scaled by reduced Bohr radius
    extent: tuple[float, float, float, float] = (
        float(np.min(Xg) / a_mu),
        float(np.max(Xg) / a_mu),
        float(np.min(Zg) / a_mu),
        float(np.max(Zg) / a_mu),
    )

    A = P
    im = ax.imshow(A, extent=extent, origin="lower", aspect="equal", cmap=cmap)

    # Axis labels
    x_z_units = r"a_\mu" if wf.use_reduced_mass else r"a_0"
    ax.set_xlabel(rf"$x / {x_z_units}$", fontsize=43, color=text_color)
    ax.set_ylabel(rf"$z / {x_z_units}$", fontsize=45, color=text_color)
    ax.xaxis.set_label_coords(x=0.5, y=-0.075)
    ax.yaxis.set_label_coords(x=-0.08, y=0.5)

    # Title and subtitle
    ax.set_title(
        "Hydrogen Wavefunction - Probability Density",
        pad=130, fontsize=44, loc="left", color=title_color
    )

    fig.text(
        x=ax.get_position().x0 + 0.07, y=0.868,
        s=r"$|\psi_{n\ell m}(r,\theta,\phi)|^{2} = |R_{n\ell}(r) Y_{\ell}^{m}(\theta,\phi)|^2$",
        fontsize=40, color=title_color
    )

    # Colormap colorbar
    cbar = plt.colorbar(im, fraction=0.046, pad=0.025)
    cbar.set_label(r"Probability density $|\psi|^{2}$ [m$^{-3}$]", fontsize=40, color=text_color, labelpad=34)
    cbar.ax.tick_params(labelsize=26, colors=text_color)
    cbar.ax.set_frame_on(not use_dark_theme)

    sf = mticker.ScalarFormatter(useMathText=True)
    sf.set_powerlimits((0, 0))
    cbar.formatter = sf
    cbar.update_ticks()

    fig.canvas.draw()
    off = cbar.ax.yaxis.get_offset_text()
    offset_str = off.get_text()
    off.set_visible(False)
    cbar.ax.text(
        0.7, 1.02, offset_str,
        transform=cbar.ax.transAxes, ha="center", va="bottom",
        fontsize=28, color=text_color
    )

    # Quantum numbers (n,l,m) label
    h, w = A.shape
    patch = A[max(h - 40, 0):h, 0:min(40, w)]  # Sample top-left corner patch to gauge brightness
    patch_val = np.nanmean(A) if patch.size == 0 else float(np.nanmean(patch))

    r, g, b, _ = im.cmap(im.norm(patch_val))  # Compute luminance of patch color
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    qn_color = "#ffffff" if luminance < 0.5 else "#000000"
    qn_outline = "#000000" if qn_color == "#ffffff" else "#ffffff"

    ax.text(
        x=0.04, y=0.95, s=f"({wf.n}, {wf.l}, {wf.m})",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=42, color=qn_color,
        path_effects=[pe.withStroke(linewidth=3.0, foreground=qn_outline)]
    )

    # Save and display figure
    ts = datetime.now().strftime("%Y%m%d%H%M%S")
    filename = f"({wf.n},{wf.l},{wf.m})_{ts}"
    plt.savefig(filename)
    plt.show()
