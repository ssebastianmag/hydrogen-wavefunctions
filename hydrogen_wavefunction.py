"""
File: hydrogen_wavefunction.py
Description: Vectorized computational core for hydrogenic bound-state eigenfunctions.

Model computations:
   - Normalized radial functions with stable log-gamma normalization & complex spherical harmonics.
   - Stationary-state wavefunction on an x–z plane grid (y=0).
   - Probability densities and radial probability distributions.
   - Reduced-mass Bohr radius and electron–nucleus reduced mass.

Model assumptions:
   - Non-relativistic, point nucleus, Schrödinger hydrogenic Hamiltonian with Coulomb potential.
   - No spin/fine-structure, external fields, or finite-nuclear-size effects.
   - Shapes broadcast; R is real-valued, Y and psi are complex.

   - SI Units implied:
     · Wavefunctions psi have units of m^{-3/2}; R has units of m^{-3/2}; Y and Z are dimensionless.
     · Cartesian coordinates (x,y,z) and radial distance r are in meters; Masses are in kilograms.

Author: Sebastian Mag
Date: August 2025
Repository: https://github.com/ssebastianmag/hydrogen-wavefunctions
"""

from typing import Optional, Literal

import numpy as np
import scipy.special as sp
from scipy.constants import physical_constants, m_e, m_p


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
            - R has units of m^{-3/2}.
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


def spherical_harmonic_Ylm(
        l: int,
        m: int,
        theta: np.ndarray,
        phi: np.ndarray
):
    """ Complex spherical harmonic Y_{l,m}(theta,phi); orthonormal on S².

        Parameters:
            l (int): Orbital angular-momentum quantum number (l ≥ 0).
            m (int): Magnetic quantum number (-l ≤ m ≤ l).
            theta (np.ndarray): Polar angle(s) in radians in the interval theta ∈ [0, π].
            phi (np.ndarray): Azimuthal angle(s) in radians in the interval phi ∈ [0, 2π).

        Returns:
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
    phi_mode: Literal["plane", "constant"] = "plane"
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

        Returns:
            Xg (np.ndarray): 2D Cartesian x-coordinate grid (y=0) in meters.
            Zg (np.ndarray): 2D Cartesian z-coordinate grid (y=0) in meters.
            psi (np.ndarray): Complex-valued coordinate-space wavefunction samples psi(x,z) on y=0.
            a_mu (float): Reduced-mass Bohr radius a_μ in meters.

        Notes:
            - psi has units of m^{-3/2}.
            - If use_reduced_mass = True, evaluate the Bohr radius with an electron–nucleus reduced mass μ,
              otherwise, use invariant electron mass m_e (μ = m_e ∴ a_μ = a₀) in the infinite–mass approximation.
            - If Z>1 and use_reduced_mass = True, M must be provided.
            - If Z=1 and M is not provided, proton mass m_p is assumed.
            - Xg, Zg and psi have the shape (grid_points, grid_points).
            - If phi_mode = "plane", phi=0 for x ≥ 0, phi = π for x < 0.
            - If phi_mode = "constant", phi ≡ phi_value across the grid.
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


def compute_radial_probability_distribution(R: np.ndarray, r: np.ndarray):
    """ Retrieve Radial probability distribution P_{n,l}(r) = r^2 * |R_{n,l}(r)|^2.

        Parameters:
            R (np.ndarray): Real-valued radial eigenfunction samples R_{n,l}(r).
            r (np.ndarray): Radial coordinate(s) in meters.

        Returns:
            P_r (np.ndarray): Real-valued P_{n,l}(r) with the same shape as r.
    """
    return (r**2) * np.abs(R) ** 2


def reduced_electron_nucleus_mass(Z: int, M: Optional[float] = None):
    """ Compute electron–nucleus reduced mass μ.

        Parameters:
            Z (int): Nuclear charge number. (Z=1 for Hydrogen, Z>1 for hydrogenic ions).
            M (float): Nuclear mass in kg.

        Returns:
            float: Two-body (electron + nucleus) system reduced mass μ in kilograms.

        Notes:
            - If Z>1, M must be provided.
            - If Z=1 and M is not provided, proton mass m_p is assumed.
    """
    if M is None:
        if Z == 1:
            M = m_p
        else:
            raise ValueError("'M' must be provided if Z>1")

    return (m_e * M) / (m_e + M)


def reduced_bohr_radius(mu: float):
    """ Compute Bohr radius evaluated with a reduced mass.

        Parameters:
            mu (float): Electron–nucleus reduced mass μ in kg.

        Returns:
            float: Reduced-mass Bohr radius a_μ in meters.
    """
    a0 = physical_constants["Bohr radius"][0]
    return a0 * (m_e / mu)
