
"""
_________________________________________________________
First Unittest: uniform phase function for solid angle
should results to a known results.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Scattering import Scatterer

npts = 200

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)

Detector = Photodiode(NA                 = 1.,             # half the sphere
                      Source             = LightSource,
                      Npts               = npts)

Detector.Fourier.Plot()

Scat = Scatterer(Diameter    = 500e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Meshes      = Detector.Meshes)

Scat.Field.Parallel = np.ones( np.shape( Scat.Field.Parallel ) )


Coupling = Detector.Coupling(Scatterer = Scat, Polarization='Parallel')


TheoVal = (2*np.pi)**2     # Theoretical value is around 6.28319**2

print('Bound angle:\n \t Theta -> {0}\n \t Phi   -> {1}\n'.format(Detector.Meshes.Theta.Boundary.Degree, Detector.Meshes.Phi.Boundary.Degree))

print('Coupling     -> {0}, \nTheoretical -> {1}'.format(Coupling, TheoVal))

if TheoVal - 0.5 < Coupling < TheoVal + 0.5:
    print('Unittest 0: Passed!')
else:
    raise Exception('Unittest 0: Failed!')









# -
