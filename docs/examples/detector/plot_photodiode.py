"""
Photodiode detector
===================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.detector import Photodiode

detector = Photodiode(
    NA=0.3,
    sampling=500,
    gamma_offset=0,
    phi_offset=0,
    polarization_filter=None
)

figure = detector.plot()

_ = figure.show()

# -
