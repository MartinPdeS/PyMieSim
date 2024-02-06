"""
Far-Fields computation
======================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
# ~~~~~~~~~~~~~~~~~~~
source = Gaussian(
    wavelength=1000e-9,
    polarization_value='right',
    polarization_type='circular',
    optical_power=1,
    NA=0.3
)

# %%
# Defining the scatterer
# ~~~~~~~~~~~~~~~~~~~~~~
scatterer = Sphere(
    diameter=1500e-9,
    source=source,
    index=1.4,
    n_medium=1.0
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
