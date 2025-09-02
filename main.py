"""
File: main.py
Description: Application entry-point to visualize hydrogenic wavefunction states.
"""

from hwf_plots import WaveFunction, plot_hydrogen_wavefunction_xz


if __name__ == "__main__":

    # (2,0,0)
    wf = WaveFunction(n=2, l=0, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=0.35)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=0.35)

    # (3,0,0)
    wf = WaveFunction(n=3, l=0, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=0.5)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=0.5)

    # (2,1,0)
    wf = WaveFunction(n=2, l=1, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2)

    # (3,1,0)
    wf = WaveFunction(n=3, l=1, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=1.8)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=1.8)

    # (3,1,1)
    wf = WaveFunction(n=3, l=1, m=1)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=1.8)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=1.8)

    # (2,1,1)
    wf = WaveFunction(n=2, l=1, m=1)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2)

    # (3,2,0)
    wf = WaveFunction(n=3, l=2, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2.2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2.2)

    # (3,2,1)
    wf = WaveFunction(n=3, l=2, m=1)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2)

    # (3,2,2)
    wf = WaveFunction(n=3, l=2, m=2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2)

    # (4,0,0)
    wf = WaveFunction(n=4, l=0, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=0.5)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=0.5)

    # (4,1,0)
    wf = WaveFunction(n=4, l=1, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=1.5)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=1.5)

    # (4,1,1)
    wf = WaveFunction(n=4, l=1, m=1)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=1.5)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=1.5)

    # (4,2,0)
    wf = WaveFunction(n=4, l=2, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=1.8)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=1.8)

    # (4,2,1)
    wf = WaveFunction(n=4, l=2, m=1)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=1.8)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=1.8)

    # (4,2,2)
    wf = WaveFunction(n=4, l=2, m=2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=1.8)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=1.8)

    # (4,3,0)
    wf = WaveFunction(n=4, l=3, m=0)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2.2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2.2)

    # (4,3,1)
    wf = WaveFunction(n=4, l=3, m=1)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2.2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2.2)

    # (4,3,2)
    wf = WaveFunction(n=4, l=3, m=2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2.2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2.2)

    # (4,3,3)
    wf = WaveFunction(n=4, l=3, m=3)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", k=2.2)
    plot_hydrogen_wavefunction_xz(wf, colormap="rocket", use_dark_theme=True, k=2.2)
