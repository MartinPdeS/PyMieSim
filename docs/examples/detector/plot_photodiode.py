"""
Photodiode detector
===================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import Photodiode

detector = Photodiode(
    NA=0.3,
    sampling=500,
    gamma_offset=-45,
    phi_offset=20,
    polarization_filter=None
)

figure = detector.plot()

_ = figure.show()

# -
