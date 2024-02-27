"""
S1 S2 function computation
==========================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.scatterer import Sphere
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
# ~~~~~~~~~~~~~~~~~~~
source = Gaussian(
    wavelength=450e-9,
    polarization_value=0,
    polarization_type='linear',
    optical_power=1,
    NA=0.3
)


# %%
# Defining the scatterer
# ~~~~~~~~~~~~~~~~~~~~~~
scatterer = Sphere(
    diameter=6e-9,
    source=source,
    index=1.4
)

# %%
# Computing the data
# ~~~~~~~~~~~~~~~~~~
data = scatterer.get_s1s2(sampling=200)

# %%
# Plotting the data
# ~~~~~~~~~~~~~~~~~
figure = data.plot()

_ = figure.show()

# -
