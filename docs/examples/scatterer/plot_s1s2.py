"""
S1 S2 function computation
==========================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.scatterer import Sphere
from PyMieSim.source import PlaneWave

# %%
# Defining the source
# ~~~~~~~~~~~~~~~~~~~
source = PlaneWave(
    wavelength=450e-9,
    polarization=0,
    amplitude=1
)

# %%
# Defining the scatterer
# ~~~~~~~~~~~~~~~~~~~~~~
scatterer = Sphere(
    diameter=600e-9,
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
