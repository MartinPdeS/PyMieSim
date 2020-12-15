
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Fields import Source

LightSource = Source(Wavelength   = 550e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 100e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Npts        = 101,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [0,180],
                 PhiOffset     = 0.0001)



Scat.Field.Plot()

Scat.S1S2.Plot()

#Scat.SPF.Plot()

Scat.Meshes.Plot()

plt.show()





# -
