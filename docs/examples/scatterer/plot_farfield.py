"""
Far-Fields computation
======================

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
# Defining the scatterer
# ~~~~~~~~~~~~~~~~~~~~~~
scatterer = Sphere(
    diameter=1500e-9,
    source=source,
    index=1.4
)

# %%
# Computing the data
# ~~~~~~~~~~~~~~~~~~
data = scatterer.get_far_field(sampling=100)

# %%
# Plotting the data
# ~~~~~~~~~~~~~~~~~
figure = data.plot()

_ = figure.show()

# -
