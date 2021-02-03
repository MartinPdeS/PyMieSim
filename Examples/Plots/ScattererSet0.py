
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Sets import ScattererSet, ExperimentalSet

LightSource = Source(Wavelength   = 950e-9,
                     Polarization = 0)



Detector0 = LPmode(NA                = 0.1,
                   Sampling          = 401,
                   GammaOffset       = 40,
                   PhiOffset         = 0,
                   Mode              = (0,1),
                   CouplingMode      = 'Mean')

Detector1 = LPmode(NA                = 0.1,
                   Sampling          = 401,
                   GammaOffset       = 40,
                   PhiOffset         = 0,
                   Mode              = (1,1),
                   CouplingMode      = 'Mean')





ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 500e-9, 100),
                       RIList        = np.linspace(1.5, 1.5, 1).round(1),
                       Source        = LightSource)





Set = ExperimentalSet(ScattererSet  = ScatSet,
                      Detectors     = [Detector0, Detector1])


Data = Set.DataFrame

Data.Plot(y='Coupling')
# -
