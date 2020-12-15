
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.functions.converts import Angle2Direct, Direct2Angle
from PyMieCoupling.classes.Detector import fiber, LPmode
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Scattering import Scatterer

PolarizationList = np.linspace(0,180,100)
CouplingLP01 = []
CouplingLP11 = []

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 40)

LP11 = LPmode(Mode          = (1, 1),
              Orientation   = 'v',
              Source        = LightSource,
              Npts          = 101,
              ThetaOffset   = 0,
              PhiOffset     = 40,
              Filter        = 'None',
              NA            = 0.2)


LP01 = LPmode(Mode          = (0, 1),
              Orientation   = 'v',
              Source        = LightSource,
              Npts          = 101,
              ThetaOffset   = 0,
              PhiOffset     = 40,
              Filter        = 0,
              NA            = 0.2)



for Polarization in PolarizationList:

    LightSource = Source(Wavelength   = 450e-9,
                         Polarization = Polarization)


    Scat = Scatterer(Diameter    = 2000e-9,
                     Source      = LightSource,
                     Index       = 1.4,
                     Meshes      = LP11.FarField.Meshes)


    CouplingLP01.append( LP01.Coupling(Scat, Mode='Centered'))
    CouplingLP11.append( LP11.Coupling(Scat, Mode='Centered'))


plt.figure(figsize=(10,5))
plt.plot(PolarizationList, CouplingLP01, 'C0', label='LP11')
plt.plot(PolarizationList, CouplingLP11, 'C1', label='LP11')
plt.grid()
plt.legend()
plt.show()
# -
