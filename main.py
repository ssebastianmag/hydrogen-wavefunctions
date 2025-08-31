"""
File: main.py
Description: Application entry-point to visualize hydrogenic wavefunction states.
"""

from hwf_plots import WaveFunction, plot_hydrogen_wavefunction_xz


if __name__ == "__main__":

    wavefunction = WaveFunction(n=3, l=2, m=1, Z=1, phi_mode="plane", extent_a_mu=20, grid_points=600)
    plot_hydrogen_wavefunction_xz(wavefunction, colormap="rocket", use_dark_theme=True)

