
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode,Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Scattering import Scatterer

PolarizationList = np.linspace(0,360,100)
Coupling = []
CouplingLP11 = []

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

LP01 = LPmode(Mode          = (0, 1,'v'),
              Source        = LightSource,
              Sampling      = 501,
              GammaOffset   = 0,
              PhiOffset     = 0,
              Filter        = 90,
              NA            = 0.2,
              CouplingMode  = 'Centered' )

Detector = Photodiode(NA                = 0.2,
                      Source            = Source,
                      Sampling          = 1001,
                      GammaOffset       = 0,
                      Filter            = 'None',
                      PhiOffset         = 85)

for Polarization in PolarizationList:

    LightSource = Source(Wavelength = 450e-9, Polarization = Polarization)


    Scat = Scatterer(Diameter    = 800e-9,
                     Source      = LightSource,
                     Index       = 1.4,
                     Meshes      = Detector.Meshes)


    Coupling.append( Detector.Coupling(Scat) )



plt.figure(figsize=(10,5))
plt.plot(PolarizationList, Coupling, 'C0', label='LP01')

#plt.xlim([0, np.max(Coupling)])
plt.grid()
plt.legend()
plt.show()
# -
