"""
LP01 Mode detector
==================

"""

# %%
# Importing the package: PyMieSim
from PyMieSim.single.detector import LPMode

detector = LPMode(
    mode_number="LP02",
    sampling=500,
    NA=0.5,
    gamma_offset=0,
    phi_offset=40,
    coupling_mode='Point'
)

figure = detector.plot()

figure.background_color = 'black'

figure.unit_size = (1200, 1200)

_ = figure.show()

# -
