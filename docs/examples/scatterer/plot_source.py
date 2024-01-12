"""
Source plottings
================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave

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
# Plotting the source
# ~~~~~~~~~~~~~~~~~~~
figure = source.plot()

_ = figure.show()

# -
