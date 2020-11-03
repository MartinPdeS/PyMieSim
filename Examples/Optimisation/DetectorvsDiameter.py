
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
from PyMieCoupling.classes.Misc import Source

npts=201

cuda = False

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Detector = Photodiode(NA                = 0.3,
                      Source            = LightSource,
                      Npts              = npts,
                      ThetaOffset       = 10,
                      PhiOffset         = 10,
                      cuda              = cuda)

Detector.Fourier.Plot('Polar')

DiameterList = np.linspace(100,3000,20).round(3) * 1e-9

DataFrame = CouplingStat(RIList       = [1.3],
                         DiameterList = DiameterList,
                         Detector     = Detector,
                         Source       = LightSource,
                         Npts         = npts,
                         cuda         = cuda)

DataFrame.plot(y='Coupling')

plt.show()








# -
