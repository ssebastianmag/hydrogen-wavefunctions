# --- --- --- --- --- --- --- --- ---
# Hydrogen Atom - Wavefunction and Electron Density Visualization

# Sebastian Mag | August 2023
# https://github.com/ssebastianmag/hydrogen-wavefunctions

# Modeling and visualization of hydrogen atom wavefunctions
# and electron probability density.
# --- --- --- --- --- --- --- --- ---

from scipy.constants import physical_constants
import matplotlib.pyplot as plt
import scipy.special as sp
import seaborn as sns
import numpy as np


# Normalized radial function Rnl(r)
def radial_function(n, l, r, a0):
    """ Compute the normalized radial part of the wavefunction using
    Laguerre polynomials and an exponential decay factor.

    Args:
        n (int): principal quantum number
        l (int): azimuthal quantum number
        r (numpy.ndarray): radial coordinate
        a0 (float): scaled Bohr radius
    Returns:
        numpy.ndarray: wavefunction radial component
    """

    # Laguerre polynomials describe how the electron density
    # changes as the distance from the nucleus increases
    laguerre = sp.genlaguerre(n - l - 1, 2 * l + 1)

    # Normalized radial distance from the nucleus
    p = 2 * r / (n * a0)

    # This factor ensures the radial wavefunction is normalized
    constant_factor = np.sqrt(
        ((2 / n * a0) ** 3 * (sp.factorial(n - l - 1))) /
        (2 * n * (sp.factorial(n + l)))
    )

    # The radial part of the wavefunction is constructed by the product of:
    # - Constant factor:
    #   Normalizes the radial wavefunction

    # - Exponential decay factor: np.exp(-p / 2)
    #   Reflects the decrease in probability of finding an
    #   electron as it moves away from the nucleus

    # - Power-law dependence on radial distance: p ** l
    #   Introduces a dependency based on the azimuthal quantum number 'l',
    #   indicating different radial behaviors for different orbitals

    # - Laguerre polynomial: laguerre(p)
    #   Captures oscillations in the electron density
    #   as a function of radial distance
    return constant_factor * np.exp(-p / 2) * (p ** l) * laguerre(p)


# Normalized angular function Ylm(θ,φ)
def angular_function(m, l, theta, phi):
    """ Compute the normalized angular part of the wavefunction using
    Legendre polynomials and a phase-shifting exponential factor.

    Args:
        m (int): magnetic quantum number
        l (int): azimuthal quantum number
        theta (numpy.ndarray): polar angle
        phi (int): azimuthal angle
    Returns:
        numpy.ndarray: wavefunction angular component
    """

    # Legendre polynomials describe the spatial arrangement and directional
    # characteristics of electron probability densities
    legendre = sp.lpmv(m, l, np.cos(theta))

    # This factor ensures that the angular wavefunction is normalized
    constant_factor = ((-1) ** m) * np.sqrt(
        ((2 * l + 1) * sp.factorial(l - np.abs(m))) /
        (4 * np.pi * sp.factorial(l + np.abs(m)))
    )

    # The angular part of the wavefunction is constructed by the product of:
    # - Constant factor:
    #   Normalizes the angular wavefunction

    # - Legendre polynomial:
    #   Describes the angular dependence of the wavefunction based on the quantum numbers.
    #   Providing insight into the orientation and shape of electron orbitals
    #   around the nucleus for given quantum numbers

    # - Exponential factor: np.real(np.exp(1.j * m * phi))
    #   Introduces a phase shift dependent on the magnetic quantum
    #   number 'm' and the azimuthal angle 'phi'
    return constant_factor * legendre * np.real(np.exp(1.j * m * phi))


# Normalized wavefunction Ψnlm(r,θ,φ) as a product of Rnl(r).Ylm(θ,φ)
def compute_wavefunction(n, l, m, a0_scale_factor):
    """ Compute the normalized wavefunction as a product
    of its radial and angular components.

    Args:
        n (int): principal quantum number
        l (int): azimuthal quantum number
        m (int): magnetic quantum number
        a0_scale_factor (float): Bohr radius scale factor
    Returns:
        numpy.ndarray: wavefunction
    """

    # The Bohr radius sets the scale of the wavefunction and determines the size of the atom.
    # By scaling it, we adapt the wavefunction's spatial extent for effective visualization
    a0 = a0_scale_factor * physical_constants['Bohr radius'][0] * 1e+12

    # Establish a grid in the z-x plane, allowing the wavefunction to assign a probability
    # value to each point. This grid aids in visualizing the electron's spatial distribution
    grid_extent = 480
    grid_resolution = 680
    z = x = np.linspace(-grid_extent, grid_extent, grid_resolution)
    z, x = np.meshgrid(z, x)

    # Using an epsilon value to prevent division by zero during the calculation of angles
    eps = np.finfo(float).eps

    # Compute the wavefunction by multiplying the radial and angular parts.
    # The radial part considers the distance from the nucleus, whereas the angular part
    # looks into the spatial orientation. Together, they define the electron's behavior
    # in the atom's vicinity
    psi = radial_function(
        n, l, np.sqrt((x ** 2 + z ** 2)), a0
    ) * angular_function(
        m, l, np.arctan(x / (z + eps)), 0
    )

    # Return the computed wavefunction, which encapsulates the quantum state
    # of an electron in a hydrogen atom. The wavefunction contains complex amplitudes
    # that provide information about the quantum state's magnitude and phase
    return psi


# Probability density |Ψ|^2
def compute_probability_density(psi):
    """ Compute the probability density of a given wavefunction.
    Args:
        psi (numpy.ndarray): wavefunction
    Returns:
        numpy.ndarray: wavefunction probability density
    """

    # Return the computed probability density, which gives the likelihood of finding
    # the electron at a specific point in space for the given quantum state. The
    # values represent the square magnitude of the wavefunction, encapsulating the
    # probability of the electron's presence in different regions of the atom
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

    # Create a new figure with specified dimensions
    fig, ax = plt.subplots(figsize=(16, 16.5))
    plt.subplots_adjust(top=0.82)
    plt.subplots_adjust(right=0.905)
    plt.subplots_adjust(left=-0.1)

    # Compute and visualize the wavefunction probability density
    # - By taking the square root of the probability density we reduce the dynamic range
    #   of the visualization, spreading out the values and making the electron's presence
    #   in low and medium probability regions more distinguishable
    psi = compute_wavefunction(n, l, m, a0_scale_factor)
    prob_density = compute_probability_density(psi)

    # Here we transpose the array to align the calculated z-x plane with Matplotlib's y-x imshow display
    im = ax.imshow(np.sqrt(prob_density).T, cmap=sns.color_palette(colormap, as_cmap=True))

    # Add a colorbar
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

    # As we examine the electron density plots corresponding to the quantum numbers above, we notice
    # that with increasing principal quantum number (n), the complexity of the wavefunction grows
    # Specifically:
    # - The number of nodes (regions where the probability density is zero) increases.
    # - The electron's spatial distribution expands, covering larger regions around the nucleus.
    # - The overall shape of the atomic orbital becomes more intricate and detailed.

    plot_wf_probability_density(9, 6, 1, 0.04, True, colormap='mako')
    plot_wf_probability_density(20, 10, 5, 0.01, True, colormap='mako')

    # For extremely high quantum numbers, the following effects can be observed:
    # - The complexity increases even further, resulting in numerous nodes and intricate patterns.
    # - Evaluating the wavefunction over a vast spatial domain becomes computationally intensive.
    # - Visualization can become cluttered, making it harder to discern specific details or features.
