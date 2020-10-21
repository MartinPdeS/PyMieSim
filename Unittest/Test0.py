
"""
_________________________________________________________
First Unittest: uniform phase function for solid angle
should results to a known results.
_________________________________________________________
"""

import numpy as np
from src.classes.Detector import Detector
from src.classes.Scattering import Scatterer
from src.functions.couplings import PointFieldCoupling

npts=201

Detector = Detector(size       = 100e-6,
                    shape      = 'circle',
                    wavelength = 400e-9,
                    npts       = npts,
                    )

Detector.magnificate( magnification=1.5 )

DiameterList = np.linspace(100,1000,50) * 1e-9

CouplingLP01, CouplingLP11 = [], []

Scat = Scatterer(diameter    = 500e-9,
                 wavelength  = 400e-9,
                 index       = 1.4,
                 npts        = npts,
                 ThetaBound  = Detector.ThetaBound,
                 ThetaOffset = 0,
                 PhiBound    = Detector.PhiBound,
                 PhiOffset   = 0
                 )

Coupling = PointFieldCoupling(Detector = Detector,
                              Source   = np.ones( np.shape( Scat.Field.Parallel ) ),  #Uniforme sphere
                              Mesh     = Scat.Meshes
                              )


TheoVal, Delta = 0.127646, 0.05  # Theoretical value is around 0.127646


print('Bound angle:\n \t Theta -> {0}\n \t Phi -> {1}\n'.format(Detector.ThetaBound, Detector.PhiBound))

print('Coupling -> {0}, \nTheoretical -> {1}'.format(Coupling, TheoVal))

if 0.127646 - Delta < Coupling < 0.127646 + Delta:
    print('Unittest 0: Passed!')
else:
    raise Exception('Unittest 0: Failed!')









# -
