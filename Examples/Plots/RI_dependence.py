
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.functions.Optimization import PlotRI
from PyMieCoupling.classes.Fields import Source



LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)


npts = 401

Detector = Photodiode(NA                = 0.08,
                      Source            = LightSource,
                      Npts              = npts,
                      ThetaOffset       = 0,
                      PhiOffset         = 50)


RIList = np.linspace(1.3, 1.6, 4).round(4)

Array = PlotRI(Diameter     = 1000e-9,
               RIList       = RIList,
               Detector     = Detector,
               Source       = LightSource,
               Npts         = npts,
               QuietMode    = True,)












# -
