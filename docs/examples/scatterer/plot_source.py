"""
Source plottings
================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.source import Gaussian

# %%
# Defining the source
# ~~~~~~~~~~~~~~~~~~~
source = Gaussian(
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
