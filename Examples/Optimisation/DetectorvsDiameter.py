
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
For different scatterer diameters.
_________________________________________________________
"""
import matplotlib.pyplot as plt
import numpy as np
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.Optimization import CouplingStat

npts=201

Detector = Photodiode(NumericalAperture = 0.3,
                      Wavelength        = 400e-9,
                      Npts              = npts,
                      ThetaOffset       = 10,
                      PhiOffset         = 10)

Detector.PlotPolar()

DiameterList = np.linspace(100,3000,20).round(3) * 1e-9

DataFrame = CouplingStat(RIList       = [1.3],
                         DiameterList = DiameterList,
                         Detector     = Detector,
                         Wavelength   = 400e-9,
                         Npts         = 101)

DataFrame.plot(y='Coupling')

plt.show()








# -
