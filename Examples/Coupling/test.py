import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Sets import ScattererSet, ExperimentalSet
from PyMieCoupling.classes.Scattering import Scatterer
from mayavi.mlab import *
from ai import cs
from mayavi import mlab

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

Detector0 = Photodiode(NA                = 0.3,
                       Sampling          = 801,
                       GammaOffset       = 0,
                       PhiOffset         = 0,
                       CouplingMode      = 'Centered')

Detector1 = LPmode(NA                = 0.1,
                   Sampling          = 401,
                   GammaOffset       = 0,
                   PhiOffset         = 0,
                   Mode              = (0,1),
                   CouplingMode      = 'Centered')

Scat = Scatterer(Index=1.4, Source = LightSource, Diameter=1000e-9)


SPF = Scat.SPF(Num=100)

SPF.Plot()
