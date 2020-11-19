
"""
_________________________________________________________
Benchmark for field coupling cython Vs numpy
_________________________________________________________
"""
import time
from PyMieCoupling.cpp.Fields import Coupling
import numpy as np


NPTS = 11


from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.Scattering import Scatterer
Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)


LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)


LP01 = LPmode(Fiber         = Fiber,
              Mode          = (0, 1),
              Source        = LightSource,
              Npts          = NPTS,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Magnification = 1)

Scat = Scatterer(Diameter    = 10e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Meshes      = LP01.Meshes)




print()

t0 = time.time() #__________________________________________________

Parallel, Perpendicular = Coupling(Scat.Index,
                                   Scat.SizeParam,
                                   LP01.Meshes.Phi.Vector.Radian.tolist(),
                                   LP01.Meshes.Theta.Vector.Radian.tolist(),
                                   LP01.Fourier.Array.tolist())
t1 = time.time() #__________________________________________________
print('Time for cython computing:============= {0}'.format(t1-t0))
print(Parallel, Perpendicular)





t0 = time.time() #__________________________________________________




res = LP01.Coupling(Scat)
t1 = time.time() #__________________________________________________
print('Time for numpy computing:============= {0}'.format(t1-t0))
print(res)



print(Scat.SizeParam)

# -
