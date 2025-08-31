"""
File: hwf_plots.py
Description: Plotting utilities for hydrogenic wavefunctions.
"""

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
        use_dark_theme: bool = False
):
    """ Plot hydrogen wavefunction restricted to the y=0 (x–z) plane.

        Parameters:
            wf (WaveFunction): Wavefunction parameters.
            colormap (str): Seaborn colormap name.
            use_dark_theme (bool): Plot theme rendering mode.
    """
    try:
        _ = sns.color_palette(colormap)
    except Exception:
        raise ValueError(f"(!) {colormap!r} is not a recognized Seaborn colormap.")

    pass