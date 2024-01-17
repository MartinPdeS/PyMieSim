"""
Photodiode detector
===================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import IntegratingSphere

detector = IntegratingSphere(
    sampling=500,
    polarization_filter=None
)

figure = detector.plot()

_ = figure.show()

# -
