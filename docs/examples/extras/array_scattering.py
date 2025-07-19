"""
Array-based scattering calculations
===================================

This example demonstrates how to compute far fields, S1/S2 amplitudes and
Stokes parameters using arbitrary ``phi`` and ``theta`` arrays.
"""

import numpy as np
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import PlaneWave
from PyMieSim.units import nanometer, degree, RIU, volt, meter

# %%
# Create a simple plane wave source
source = PlaneWave(
    wavelength=632.8 * nanometer,
    polarization=0 * degree,
    ampliutde=1 * volt / meter,
)

# %%
# Define the scatterer
scatterer = Sphere(
    diameter=200 * nanometer,
    property=1.5 * RIU,
    medium_property=1.0 * RIU,
    source=source,
)

# Define arbitrary angle arrays
phi = np.linspace(0, np.pi, 8)
theta = np.linspace(0, np.pi / 2, 8)

# %%
# Far-field complex fields
E_para, E_perp = scatterer.get_far_field_array(phi, theta)
print(E_para.shape, E_perp.shape)

# %%
# S1 and S2 scattering amplitudes
S1, S2 = scatterer.get_s1s2_array(phi)
print(S1.shape, S2.shape)

# %%
# Stokes parameters
I, Q, U, V = scatterer.get_stokes_array(phi, theta)
print(I.shape, Q.shape, U.shape, V.shape)
