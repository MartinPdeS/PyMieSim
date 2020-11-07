
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""
import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.functions.Optimization import ComputeSTD

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

npts=51

cuda = True

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)


Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = LPmode(Fiber         = Fiber,
              Mode          = (1, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Magnification = 1,
              cuda          = cuda)

LP01 = LPmode(Fiber         = Fiber,
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = npts,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Magnification = 1,
              cuda          = cuda)

STD = ComputeSTD(RIList       = np.linspace(1.3, 1.8, 20).tolist(),
                 DiameterList = np.linspace(1e-6, 5e-6, 10).tolist(),
                 Detector     = LP01,
                 QuietMode    = False,
                 cuda         = cuda,
                 Source       = LightSource)


print(STD)





# -
