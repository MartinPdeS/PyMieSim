"""
Scatterer footprint
===================

"""

# %%
# Importing the package: PyOptik, PyMieSim
from PyMieSim.scatterer import Sphere
from PyMieSim.detector import LPmode
from PyMieSim.source import PlaneWave
from PyOptik import ExpData

# %%
# Defining the source
# ~~~~~~~~~~~~~~~~~~~
source = PlaneWave(
    wavelength=1000e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1,
    NA=0.3
)

# %%
# Defining the scatterer
# ~~~~~~~~~~~~~~~~~~~~~~
scatterer = Sphere(
    diameter=2000e-9,
    source=source,
    material=ExpData('BK7')
)

# %%
# Defining the detector
# ~~~~~~~~~~~~~~~~~~~~~
detector = LPmode(
    mode_number="LP21",
    NA=0.3,
    sampling=200,
    gamma_offset=0,
    phi_offset=0,
    coupling_mode='Point'
)

# %%
# Computing the data
# ~~~~~~~~~~~~~~~~~~
data = detector.get_footprint(scatterer)

# %%
# Plotting the data
# ~~~~~~~~~~~~~~~~~
figure = data.plot()

_ = figure.show()

# -
