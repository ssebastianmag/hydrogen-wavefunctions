# Hydrogen Wavefunctions - Probability density plots

Modeling and visualization of solutions 
for the hydrogen atom wavefunction and 
electron probability density plots.

* Python 3.10.4
* Matplotlib 3.6.0
* Seaborn 0.12.0
* NumPy 1.23.3
* SciPy 1.9.1

---
### Files / Modules - Execution alternatives
* [IPython Notebook / Jupyter Notebook: hydrogen_wavefunctions.ipynb](hydrogen_wavefunctions.ipynb)
* [Module with command line arguments: hydrogen_wavefunctions_run.py](hydrogen_wavefunctions_run.py)
* [Main module: hydrogen_wavefunctions.py](hydrogen_wavefunctions.py)
---

## Content
#### Theory
* [QM Introduction: wavefunctions, 
atomic orbitals and probability](#qm-introduction-wavefunctions-atomic-orbitals-and-probability)
* [Quantum numbers](#quantum-numbers)
* [Hydrogen wavefunction modeling](#hydrogen-wavefunction-modeling)
#### Implementation
* [Execution and examples](#execution)
---

## Theory
### QM Introduction: wavefunctions, atomic orbitals and probability

Quantum mechanics use wavefunctions to describe the mathematical 
relationship between the motion of electrons in atoms and molecules 
and their energies.

A wavefunction (Ψ) is a mathematical function that relates the 
location of an electron at a given point in space 
(identified by x, y, and z coordinates) to the amplitude of its 
wave, which corresponds to its energy.

<br>
<p align='center'>
  <img src='img/Hydrogen Wavefunction Probability density plots.png' width=90% />
</p>
<p align='center'>
    <i>Electron probability density in orbital cross-sections</i>
</p>

The probability of finding an electron at a given point is proportional to 
the square of the wavefunction at that point, leading to a distribution of 
probabilities in space that is often portrayed as an electron density plot.

An atomic orbital is a mathematical function describing 
the wave-like behavior of an electron in an atom. Orbitals are 
mathematically derived regions of space where the electron can be calculated to be present.

<br>
<p align='center'>
    <img src='img/wavefunction (4,3,1)-dp.png' width=50% />
</p>
<p align='center'>
    <i>Light shaded areas represent high probability density, 
    darker areas represent lower probability density</i>
</p>

The description of electron distribution as standing waves leads naturally
to the existence of sets of quantum numbers (n, l, m) characteristic of each 
wavefunction.

---

### Quantum numbers

Schrödinger’s approach uses three quantum numbers (n, l, m) 
to specify any wavefunction. The quantum numbers provide information 
about the spatial distribution of an electron. 

Although n can be any positive integer, only certain values of 
l and m are allowed for a given value of n.

* **Principal quantum number (n):** `( 1 <= n )`

Average relative distance of 
an electron from the nucleus.

* **Azimuthal quantum number (l):** `( 0 <= l <= n-1 )`

Shape of the region of 
space occupied by the electron.

* **Magnetic quantum number (m):** `( -l <= m <= l )`

Orientation of the region in space occupied by an electron with respect to 
an applied magnetic field.

---

### Hydrogen wavefunction modeling

We may solve Schrödinger’s equation more easily if we express it in 
terms of the spherical coordinates (r, θ, φ) instead of rectangular 
coordinates (x, y, z).

In spherical coordinates, the variable r is the radial coordinate, 
θ is the polar angle (relative to the vertical z-axis), 
and φ is the azimuthal angle (relative to the x-axis).

<p align='center'>
  <img src='img/coordinate_system.png' width=38% />
</p>
<p align='center'>
    <i>Relationship between the spherical and rectangular coordinate systems</i>
</p>

The wavefunctions for the hydrogen atom depend upon the three variables r, θ, and φ 
and the three quantum numbers n, l, and m. 
The solutions to the hydrogen atom Schrödinger equation are functions that are 
products of a Radial Function and a Spherical Harmonic Function.

<br>
<p align='center'>
  <img src='img/normalized_wf.png' width=60%/>
</p>
<p align='center'>
    <i>Normalized position wavefunctions, given in spherical coordinates</i>
</p>

<p align='center'>
    <img src='img/radial_function.png' width=60%/>
</p>
<p align='center'>
    <i>Normalized radial function</i>
</p>

<p align='center'>
    <img src='img/angular_function.png' width=60%/>
</p>
<p align='center'>
    <i>Normalized angular function (spherical harmonic function)</i>
</p>

The absolute square of the wavefunction, evaluated at r, θ,
and φ gives the probability density of finding the electron.

---
## Implementation

### Execution
* [IPython Notebook / Jupyter Notebook](hydrogen_wavefunctions.ipynb)
* [Module with command line arguments](hydrogen_wavefunctions_run.py)
<br><br>

#### Command line arguments:

```
$ python hydrogen_wavefunctions_run.py --help   
```

```
usage: hydrogen_wavefunctions_run.py [-h] [-dp] n l m a0      

Hydrogen Wavefunction electron probability-density plot by definition of 
quantum numbers n, l, m and a bohr radius augmentation coefficient

positional arguments:
  n           (n) Principal quantum number (dtype: integer) (1 <= n)
  l           (l) Azimuthal quantum number (dtype: integer) (0 <= l <= n-1)
  m           (m) Magnetic quantum number (dtype: integer) (-l <= m <= l)
  a0          (a0) Bohr radius augmentation coefficient (dtype: float) (0 < a0 <= 1)

options:
  -h, --help  show this help message and exit
  -dp         Dark palette: Plot with a less bright color scheme
```
---

#### Input args:
    $ python hydrogen_wavefunctions_run.py 2 1 1 0.6

|     |                       Argument                       | Value |  Constraint   |
|:---:|:----------------------------------------------------:|:-----:|:-------------:|
|  n  |             Principal quantum number (n)             |   2   |    1 <= n     |
|  l  |            Azimuthal quantum number  (l)             |   1   | 0 <= l <= n-1 |
|  m  |             Magnetic quantum number (m)              |   1   | -l <= m <= l  |
| a0  | Bohr radius augmentation coefficient (a<sub>0</sub>) |  0.6  |  0 < a0 <= 1  |

#### Output:

<p align='center'>
  <img src='img/wavefunction (2,1,1).png' width=60% />
</p>

---

#### Input args:
    $ python hydrogen_wavefunctions_run.py 2 1 1 0.6 -dp

|     |                       Argument                       | Value |  Constraint   |
|:---:|:----------------------------------------------------:|:-----:|:-------------:|
|  n  |             Principal quantum number (n)             |   2   |    1 <= n     |
|  l  |            Azimuthal quantum number  (l)             |   1   | 0 <= l <= n-1 |
|  m  |             Magnetic quantum number (m)              |   1   | -l <= m <= l  |
| a0  | Bohr radius augmentation coefficient (a<sub>0</sub>) |  0.6  |  0 < a0 <= 1  |
| -dp |   Dark palette; enable a less bright color scheme    | True  |               |

#### Output:

<p align='center'>
  <img src='img/wavefunction (2,1,1)-dp.png' width=60% />
</p>

---

#### Input args:
    $ python hydrogen_wavefunctions_run.py 3 2 1 0.3

|     |                       Argument                       | Value |  Constraint   |
|:---:|:----------------------------------------------------:|:-----:|:-------------:|
|  n  |             Principal quantum number (n)             |   3   |    1 <= n     |
|  l  |            Azimuthal quantum number  (l)             |   2   | 0 <= l <= n-1 |
|  m  |             Magnetic quantum number (m)              |   1   | -l <= m <= l  |
| a0  | Bohr radius augmentation coefficient (a<sub>0</sub>) |  0.3  |  0 < a0 <= 1  |

#### Output:

<p align='center'>
  <img src='img/wavefunction (3,2,1).png' width=60% />
</p>

---

#### Input args:
    $ python hydrogen_wavefunctions_run.py 3 2 1 0.3 -dp

|     |                       Argument                       | Value |  Constraint   |
|:---:|:----------------------------------------------------:|:-----:|:-------------:|
|  n  |             Principal quantum number (n)             |   3   |    1 <= n     |
|  l  |            Azimuthal quantum number  (l)             |   2   | 0 <= l <= n-1 |
|  m  |             Magnetic quantum number (m)              |   1   | -l <= m <= l  |
| a0  | Bohr radius augmentation coefficient (a<sub>0</sub>) |  0.3  |  0 < a0 <= 1  |
| -dp |   Dark palette; enable a less bright color scheme    | True  |               |

#### Output:

<p align='center'>
  <img src='img/wavefunction (3,2,1)-dp.png' width=60% />
</p>

---

#### Input args:
    $ python hydrogen_wavefunctions_run.py 4 3 0 0.2

|     |                       Argument                       | Value |  Constraint   |
|:---:|:----------------------------------------------------:|:-----:|:-------------:|
|  n  |             Principal quantum number (n)             |   4   |    1 <= n     |
|  l  |            Azimuthal quantum number  (l)             |   3   | 0 <= l <= n-1 |
|  m  |             Magnetic quantum number (m)              |   0   | -l <= m <= l  |
| a0  | Bohr radius augmentation coefficient (a<sub>0</sub>) |  0.2  |  0 < a0 <= 1  |

#### Output:

<p align='center'>
  <img src='img/wavefunction (4,3,0).png' width=60% />
</p>

---

#### Input args:
    $ python hydrogen_wavefunctions_run.py 4 3 0 0.2 -dp

|     |                       Argument                       | Value |  Constraint   |
|:---:|:----------------------------------------------------:|:-----:|:-------------:|
|  n  |             Principal quantum number (n)             |   4   |    1 <= n     |
|  l  |            Azimuthal quantum number  (l)             |   3   | 0 <= l <= n-1 |
|  m  |             Magnetic quantum number (m)              |   0   | -l <= m <= l  |
| a0  | Bohr radius augmentation coefficient (a<sub>0</sub>) |  0.2  |  0 < a0 <= 1  |
| -dp |   Dark palette; enable a less bright color scheme    | True  |               |

#### Output:

<p align='center'>
  <img src='img/wavefunction (4,3,0)-dp.png' width=60% />
</p>

---
