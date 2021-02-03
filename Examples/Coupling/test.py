import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Sets import ScattererSet



LightSource = Source(Wavelength   = 930e-9,
                     Polarization = 0)


Detector0 = LPmode(NA                = 0.2,
                   Sampling          = 301,
                   GammaOffset       = 0,
                   PhiOffset         = 180,
                   Mode              = (0,1))



Detector1 = LPmode(NA                = 0.2,
                   Sampling          = 301,
                   GammaOffset       = 0,
                   PhiOffset         = 180,
                   Mode              = (1,1))

DiameterList  = np.linspace(10e-9, 10000e-9, 500)
RIList        = np.linspace(1.36, 1.56, 1).round(1)

Set0 = ScattererSet(DiameterList  = DiameterList,
                    RIList        = RIList,
                    Source        = LightSource)


a = Set0.Qsca()
a.Plot(logy=True)
