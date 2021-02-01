
"""
_________________________________________________________
First Unittest: uniform phase function for solid angle
should results to a known results.
_________________________________________________________
"""

import numpy as np
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.utils import Source
npts=201

LightSource = Source(Wavelength = 400e-9, Polarization=0)


LP01 = LPmode(NA          = 0.2,
              Mode        = (0, 1),
              Sampling    = 200,
              GammaOffset = 0,
              PhiOffset   = 0,)

LP11 = LPmode(NA          = 0.2,
              Mode        = (1, 1),
              Sampling    = 200,
              GammaOffset = 0,
              PhiOffset   = 0,)


DiameterList = np.linspace(100,1000,50) * 1e-9

Scat = Scatterer(Diameter    = 500e-9,
                 Source      = LightSource,
                 Index       = 1.4)

Scat.Parallel = np.ones( np.shape( Scat.Field.Parallel ) )   #Uniforme sphere

LP11Coupling, _ = PointFieldCoupling(Detector = LP11,
                                     Source   = Scat)

LP01Coupling, _ = PointFieldCoupling(Detector = LP01,
                                     Source   = Scat)

LP01TheoVal, Delta = 1.0, 0.05   # Theoretical value is 0
LP11TheoVal, Delta = 0.0, 0.05   # Theoretical value is 0

print('Bound angle:\n \t Theta -> {0}\n \t Phi   -> {1}\n'.format\
      (LP11.Meshes.Theta.Boundary.Degree, LP11.Meshes.Phi.Boundary.Degree))

print('LP01 Coupling -> {0}, \nTheoretical -> {1}'.format(LP01Coupling, LP01TheoVal))
print('LP11 Coupling -> {0}, \nTheoretical -> {1}'.format(LP11Coupling, LP11TheoVal))

if 0.0 - Delta < LP11Coupling < 0.0 + Delta:
    print('Unittest 0: Passed!')
else:
    raise Exception('Unittest 0: Failed!')









# -
