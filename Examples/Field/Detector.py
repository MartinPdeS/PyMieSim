
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face detector
_________________________________________________________
"""

from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Detector import Photodiode
import matplotlib.pyplot as plt

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Detector = Photodiode(NA                = 0.2,
                      Source            = Source,
                      Npts              = 201,
                      ThetaOffset       = 0,
                      PhiOffset         = 0)


Detector.PhiOffset = 10

Detector.ThetaOffset = 60

Detector.FarField.Plot()

plt.show()




# -
