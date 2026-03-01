"""
Sphere: Qsca vs wavelength
==========================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy as np
from PyMieSim.units import ureg

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian, PolarizationSet
from PyMieSim.experiment import Setup
from PyOptik import MaterialBank
