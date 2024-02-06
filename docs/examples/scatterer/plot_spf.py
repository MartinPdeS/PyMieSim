"""
SPF computation
===============

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
# ~~~~~~~~~~~~~~~~~~~
source = Gaussian(
    wavelength=500e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1,
    NA=0.3
)

# %%
# Defining the scatterer
# ~~~~~~~~~~~~~~~~~~~~~~
scatterer = Sphere(
    diameter=1200e-9,
    source=source,
    index=1.0,
    n_medium=1.2,
    # index=1.2,
    # n_medium=1.0
)

# %%
# Computing the data
# ~~~~~~~~~~~~~~~~~~
data = scatterer.get_spf(sampling=300)

# %%
# Plotting the data
# ~~~~~~~~~~~~~~~~~
figure = data.plot()

_ = figure.show()

# -
