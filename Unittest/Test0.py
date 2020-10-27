
"""
_________________________________________________________
First Unittest: uniform phase function for solid angle
should results to a known results.
_________________________________________________________
"""

import numpy as np
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling

npts = 501

Detector = Photodiode(NumericalAperture  = 1.,         #half the sphere
                      Wavelength         = 400e-9,
                      Npts               = npts)

DiameterList = np.linspace(100,1000,50).round(4) * 1e-9

Scat = Scatterer(Diameter    = 500e-9,
                 Wavelength  = 1200e-9,
                 Index       = 1.4,
                 Meshes      = Detector.Meshes)

Scat.Field.Parallel = np.ones( np.shape( Scat.Field.Parallel ) )


Para, Perp = PointFieldCoupling(Detector = Detector,
                                Source   = Scat)

TheoVal, Delta = (2*np.pi)**2, 0.5  # Theoretical value is around 0.127646

print('Bound angle:\n \t Theta -> {0}\n \t Phi   -> {1}\n'.format(Detector.Meshes.Theta.Boundary.Degree, Detector.Meshes.Phi.Boundary.Degree))

print('Coupling     -> {0}, \nTheoretical -> {1}'.format(Para, TheoVal))

if TheoVal - Delta < Para < TheoVal + Delta:
    print('Unittest 0: Passed!')
else:
    raise Exception('Unittest 0: Failed!')









# -
