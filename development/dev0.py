import numpy as np
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import PlaneWave
from PyMieSim.units import nanometer, degree, watt, RIU, meter, volt

# %%
# Create a simple plane wave source
source = PlaneWave(
    wavelength=632.8 * nanometer,
    polarization=0 * degree,
    amplitude=1 * volt / meter,
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

# %%
# S1 and S2 scattering amplitudes
S1, S2 = scatterer.get_s1s2_array(phi)

print(S1)