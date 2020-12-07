
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Fields import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = None)

Scat = Scatterer(Diameter    = 60e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Npts        = 101,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [-180, 180])


Scat.S1S2.Plot()

Scat.SPF.Plot()

Scat.Stokes.Plot()

plt.show()

#print(Scat.Stokes)




# -
