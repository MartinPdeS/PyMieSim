
"""
_________________________________________________________
First Unittest: uniform phase function for solid angle
should results to a known results.
_________________________________________________________
"""

import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling

npts=201

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = mode(fiber      = Fiber,
            LPmode     = (2, 1),
            wavelength = 400e-9,
            npts       = npts,
            )


DiameterList = np.linspace(100,1000,50) * 1e-9

CouplingLP01, CouplingLP11 = [], []

Scat = Scatterer(diameter    = 500e-9,
                 wavelength  = 400e-9,
                 index       = 1.4,
                 npts        = npts,
                 Meshes      = LP11.Meshes
                 )
Scat.Field.Parallel = np.ones( np.shape( Scat.Field.Parallel ) )

Coupling = PointFieldCoupling(Detector = LP11,
                              Source   = Scat,  #Uniforme sphere
                              Field = 'Parallel'
                              )


TheoVal, Delta = 0.0, 0.05   # Theoretical value is 0

print('Bound angle:\n \t Theta -> {0}\n \t Phi -> {1}\n'.format\
      (LP11.Meshes.Theta.Boundary.Degree, LP11.Meshes.Phi.Boundary.Degree))

print('Coupling -> {0}, \nTheoretical -> {1}'.format(Coupling, TheoVal))

if 0.127646 - Delta < Coupling < 0.127646 + Delta:
    print('Unittest 0: Passed!')
else:
    raise Exception('Unittest 0: Failed!')









# -
