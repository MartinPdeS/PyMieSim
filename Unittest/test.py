import numpy as np
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Sets import ScattererSet

LightSource = Source(Wavelength   = 950e-9,
                     Polarization = 0)


ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 3500e-9, 100),
                       RIList        = np.linspace(1.5, 1.5, 1).round(1),
                       Source        = LightSource)


Qsca = ScatSet.Qsca()

Qsca.Plot()
