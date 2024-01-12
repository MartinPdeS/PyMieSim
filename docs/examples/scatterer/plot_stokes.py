"""
Stokes parameters computation
=============================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave

# %%
# Defining the source
# ~~~~~~~~~~~~~~~~~~~
source = PlaneWave(
    wavelength=750e-9,
    polarization_value='right',
    polarization_type='circular',
    optical_power=1,
    NA=0.3
)


# %%
# Defining the scatterer
# ~~~~~~~~~~~~~~~~~~~~~~
scatterer = Sphere(
    diameter=300e-9,
    source=source,
    index=1.4
)

# %%
# Computing the data
# ~~~~~~~~~~~~~~~~~~
data = scatterer.get_stokes(sampling=100)

# %%
# Plotting the data
# ~~~~~~~~~~~~~~~~~
figure = data.plot()

_ = figure.show()

# -
