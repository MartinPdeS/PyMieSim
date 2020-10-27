
"""
_________________________________________________________
First Unittest: uniform phase function for solid angle
should results to a known results.
_________________________________________________________
"""

import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling

npts=201

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP01 = LPmode(Fiber      = Fiber,
              Mode     = (0, 1),
              Wavelength = 400e-9,
              Npts       = npts,
              ThetaOffset = 0,
              PhiOffset   = 0,)

LP11 = LPmode(Fiber      = Fiber,
              Mode     = (1, 1),
              Wavelength = 400e-9,
              Npts       = npts,
              ThetaOffset = 0,
              PhiOffset   = 0,)


DiameterList = np.linspace(100,1000,50) * 1e-9

Scat = Scatterer(Diameter    = 500e-9,
                 Wavelength  = 400e-9,
                 Index       = 1.4,
                 Meshes      = LP11.Meshes
                 )

Scat.Field.Parallel = np.ones( np.shape( Scat.Field.Parallel ) )   #Uniforme sphere

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
