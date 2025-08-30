from typing import Optional, Literal

import numpy as np
import scipy.special as sp
from scipy.constants import physical_constants, m_e, m_p
from scipy import integrate as sci_integrate


def radial_wavefunction_Rnl(
        n: int,
        l: int,
        r: np.ndarray,
        Z: int = 1,
        use_reduced_mass: bool = True,
        M: Optional[float] = None
):
    """ Normalized hydrogenic radial wavefunction R_{n,l}(r).

        Parameters:
            n (int): Principal quantum number (n ≥ 1).
            l (int): Orbital angular-momentum quantum number (0 ≤ l ≤ n-1).
            r (np.ndarray): Radial coordinate(s) in meters.
            Z (int): Nuclear charge number. (Z=1 for Hydrogen, Z>1 for hydrogenic ions).
            use_reduced_mass (bool): Reduced-mass μ correction in Bohr radius.
            M (float): Nuclear mass in kg.

        Returns:
            R (np.ndarray): Real-valued radial eigenfunction samples with broadcasted shape of (r).

        Notes:
            - If use_reduced_mass = True, evaluate the Bohr radius with an electron–nucleus reduced mass μ,
              otherwise, use invariant electron mass m_e (μ = m_e ∴ a_μ = a₀) in the infinite–mass approximation.
            - If Z>1 and use_reduced_mass = True, M must be provided.
            - If Z=1 and M is not provided, proton mass m_p is assumed.
    """
    if not (n >= 1 and 0 <= l <= n - 1):
        raise ValueError("(!) Quantum numbers (n,l) must satisfy n ≥ 1 and 0 ≤ l ≤ n-1")

    mu = reduced_electron_nucleus_mass(Z, M) if use_reduced_mass else m_e
    a_mu = reduced_bohr_radius(mu)

    rho = 2.0 * Z * r / (n * a_mu)
    L = sp.eval_genlaguerre(n - l - 1, 2 * l + 1, rho)

    # Stable normalization prefactor using log-gamma
    # pref = (2Z/(n a_mu))^(3/2) * sqrt( (n-l-1)! / (2n * (n+l)!))

    log_pref = 1.5 * np.log(2.0 * Z / (n * a_mu))
    log_pref += 0.5 * (sp.gammaln(n - l) - (np.log(2.0 * n) + sp.gammaln(n + l + 1)))
    pref = np.exp(log_pref)
    R = pref * np.exp(-rho / 2.0) * np.power(rho, l) * L
    return R


def spherical_harmonic_Ylm(l: int, m: int, theta: np.ndarray, phi: np.ndarray):
    """ Complex spherical harmonic Y_{l,m}(theta,phi); orthonormal on S².

        Parameters:
            l (int): Orbital angular-momentum quantum number (l ≥ 0).
            m (int): Magnetic quantum number (-l ≤ m ≤ l).
            theta (np.ndarray): Polar angle(s) in radians in the interval theta ∈ [0, π].
            phi (np.ndarray): Azimuthal angle(s) in radians in the interval phi ∈ [0, 2π).

        Returns:s
            Y (np.ndarray): Complex-valued spherical harmonic with broadcasted shape of (theta, phi).
    """
    if not (l >= 0 and -l <= m <= l):
        raise ValueError("(!) Quantum numbers (l,m) must satisfy l ≥ 0 and -l ≤ m ≤ l")

    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    Y = sp.sph_harm_y(l, m, theta, phi)
    return Y


def compute_psi_xz_slice(
    n: int,
    l: int,
    m: int,
    Z: int = 1,
    use_reduced_mass: bool = True,
    M: Optional[float] = None,
    extent_a_mu: float = 20.0,
    grid_points: int = 600,
    phi_value: float = 0.0,
    phi_mode: Literal["plane", "constant"] = "plane",
):
    """ Evaluate psi_{n,l,m}(x,0,z), hydrogenic eigenfunction restricted to the y=0 (x–z) plane.

        Parameters:
            n (int): Principal quantum number (n ≥ 1).
            l (int): Orbital angular-momentum quantum number (0 ≤ l ≤ n-1).
            m (int): Magnetic quantum number (-l ≤ m ≤ l).
            Z (int): Nuclear charge number. (Z=1 for Hydrogen, Z>1 for hydrogenic ions).
            use_reduced_mass (bool): Reduced-mass μ correction in Bohr radius.
            M (float): Nuclear mass in kg.
            extent_a_mu (float): Half-width of the square grid in units of the reduced-mass Bohr radius.
            grid_points (int): Number of points per Cartesian axis.
            phi_value (float): Azimuth phi (radians); Only used if phi_mode="constant".
            phi_mode (str): Azimuthal prescription on y=0.

        Notes:
            - If use_reduced_mass = True, evaluate the Bohr radius with an electron–nucleus reduced mass μ,
              otherwise, use invariant electron mass m_e (μ = m_e ∴ a_μ = a₀) in the infinite–mass approximation.
            - If Z>1 and use_reduced_mass = True, M must be provided.
            - If Z=1 and M is not provided, proton mass m_p is assumed.
            - Xg, Zg (meters) and psi have the shape (grid_points, grid_points).
            - If phi_mode = "plane", phi=0 for x ≥ 0, phi = π for x < 0.
            - If phi_mode = "constant", phi ≡ phi_value across the grid.

        Returns:
            Xg (np.ndarray): 2D Cartesian x-coordinate grid (y=0).
            Zg (np.ndarray): 2D Cartesian z-coordinate grid (y=0).
            psi (np.ndarray): Complex-valued coordinate-space wavefunction samples psi(x,z) on y=0.
            a_mu (float): Reduced-mass Bohr radius a_μ in meters.
    """
    if not (n >= 1 and 0 <= l <= n - 1 and -l <= m <= l):
        raise ValueError("(!) Quantum numbers (n,l,m) must satisfy n ≥ 1, 0 ≤ l ≤ n-1, and -l ≤ m ≤ l")

    mu = reduced_electron_nucleus_mass(Z, M) if use_reduced_mass else m_e
    a_mu = reduced_bohr_radius(mu)

    r_max = extent_a_mu * a_mu
    axis = np.linspace(-r_max, r_max, grid_points)
    Zg, Xg = np.meshgrid(axis, axis, indexing="ij")

    r = np.hypot(Xg, Zg)
    cos_theta = np.empty_like(r)
    np.divide(Zg, r, out=cos_theta, where=(r > 0))
    cos_theta[r == 0] = 1.0

    theta = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    phi = np.where(Xg >= 0.0, 0.0, np.pi) if phi_mode == "plane" else np.full_like(r, float(phi_value))

    # Compute psi_{n,l,m}(r,theta,phi) = R_{n,l}(r) · Y_{l,m}(theta,phi)
    R = radial_wavefunction_Rnl(n, l, r, Z=Z, use_reduced_mass=use_reduced_mass, M=M)
    Y = spherical_harmonic_Ylm(l, m, theta, phi)
    psi = R * Y

    return Xg, Zg, psi, a_mu


def compute_probability_density(psi: np.ndarray):
    """ Retrieve probability density |psi|^2.

        Parameters:
            psi (np.ndarray): Complex-valued coordinate-space wavefunction samples psi(x,z) on y=0.

        Returns:
            P (np.ndarray): Real-valued |psi|^2 with the same shape as psi.
    """
    return np.abs(psi) ** 2


def plot_wf_probability_density(n, l, m, a0_scale_factor, dark_theme=False, colormap='rocket'):
    """ Plot the probability density of the hydrogen
    atom's wavefunction for a given quantum state (n,l,m).

    Args:
        n (int): principal quantum number, determines the energy level and size of the orbital
        l (int): azimuthal quantum number, defines the shape of the orbital
        m (int): magnetic quantum number, defines the orientation of the orbital
        a0_scale_factor (float): Bohr radius scale factor
        dark_theme (bool): If True, uses a dark background for the plot, defaults to False
        colormap (str): Seaborn plot colormap, defaults to 'rocket'
    """

    # Quantum numbers validation
    if not isinstance(n, int) or n < 1:
        raise ValueError('n should be an integer satisfying the condition: n >= 1')
    if not isinstance(l, int) or not (0 <= l < n):
        raise ValueError('l should be an integer satisfying the condition: 0 <= l < n')
    if not isinstance(m, int) or not (-l <= m <= l):
        raise ValueError('m should be an integer satisfying the condition: -l <= m <= l')

    # Colormap validation
    try:
        sns.color_palette(colormap)
    except ValueError:
        raise ValueError(f'{colormap} is not a recognized Seaborn colormap.')

    # Configure plot aesthetics using matplotlib rcParams settings
    plt.rcParams['font.family'] = 'STIXGeneral'
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['xtick.major.width'] = 4
    plt.rcParams['ytick.major.width'] = 4
    plt.rcParams['xtick.major.size'] = 15
    plt.rcParams['ytick.major.size'] = 15
    plt.rcParams['xtick.labelsize'] = 30
    plt.rcParams['ytick.labelsize'] = 30
    plt.rcParams['axes.linewidth'] = 4

    fig, ax = plt.subplots(figsize=(16, 16.5))
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(right=0.905)
    plt.subplots_adjust(left=-0.1)

    # Compute and visualize the wavefunction probability density
    psi = compute_wavefunction(n, l, m, a0_scale_factor)
    prob_density = compute_probability_density(psi)

    # Here we transpose the array to align the calculated z-x plane with Matplotlib's y-x imshow display
    im = ax.imshow(np.sqrt(prob_density).T, cmap=sns.color_palette(colormap, as_cmap=True))

    cbar = plt.colorbar(im, fraction=0.046, pad=0.03)
    cbar.set_ticks([])

    # Apply dark theme parameters
    if dark_theme:
        theme = 'dt'
        background_color = sorted(
            sns.color_palette(colormap, n_colors=100),
            key=lambda color: 0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2]
        )[0]
        plt.rcParams['text.color'] = '#dfdfdf'
        title_color = '#dfdfdf'
        fig.patch.set_facecolor(background_color)
        cbar.outline.set_visible(False)
        ax.tick_params(axis='x', colors='#c4c4c4')
        ax.tick_params(axis='y', colors='#c4c4c4')
        for spine in ax.spines.values():
            spine.set_color('#c4c4c4')

    else:  # Apply light theme parameters
        theme = 'lt'
        plt.rcParams['text.color'] = '#000000'
        title_color = '#000000'
        ax.tick_params(axis='x', colors='#000000')
        ax.tick_params(axis='y', colors='#000000')

    ax.set_title('Hydrogen Atom - Wavefunction Electron Density', pad=130, fontsize=44, loc='left', color=title_color)
    ax.text(0, 722, (
        r'$|\psi_{n \ell m}(r, \theta, \varphi)|^{2} ='
        r' |R_{n\ell}(r) Y_{\ell}^{m}(\theta, \varphi)|^2$'
    ), fontsize=36)
    ax.text(30, 615, r'$({0}, {1}, {2})$'.format(n, l, m), color='#dfdfdf', fontsize=42)
    ax.text(770, 140, 'Electron probability distribution', rotation='vertical', fontsize=40)
    ax.text(705, 700, 'Higher\nprobability', fontsize=24)
    ax.text(705, -60, 'Lower\nprobability', fontsize=24)
    ax.text(775, 590, '+', fontsize=34)
    ax.text(769, 82, '−', fontsize=34, rotation='vertical')
    ax.invert_yaxis()

    # Save and display the plot
    plt.savefig(f'({n},{l},{m})[{theme}].png')
    plt.show()


# - - - Example probability densities for various quantum states (n,l,m)
if __name__ == '__main__':

    plot_wf_probability_density(2, 1, 1, 0.6, True)
    plot_wf_probability_density(2, 1, 1, 0.6)

    plot_wf_probability_density(3, 2, 1, 0.3, True)
    plot_wf_probability_density(3, 2, 1, 0.3)

    plot_wf_probability_density(4, 3, 0, 0.2, dark_theme=True, colormap='magma')
    plot_wf_probability_density(4, 3, 0, 0.2, colormap='magma')
    plot_wf_probability_density(4, 3, 1, 0.2, dark_theme=True, colormap='mako')
