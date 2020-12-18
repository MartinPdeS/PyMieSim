
"""
_________________________________________________________
Plotting of Stokes parameter for scattere far-field.
_________________________________________________________
"""

import matplotlib.pyplot as plt
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Fields import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Scat = Scatterer(Diameter    = 100e-9,
                 Source      = LightSource,
                 Index       = 1.4,
                 Npts        = 101,
                 ThetaBound  = [-180, 180],
                 PhiBound    = [0,180],
                 PhiOffset     = 0.)




Scat.FarField.Plot()

#Scat.S1S2.Plot()

Scat.SPF.Plot()

#Scat.Meshes.Plot()


plt.show()





# -
