from PyMieSim.Source import PlaneWave
from PyMieSim.Detector import Photodiode


Source = PlaneWave(Wavelength   = 450e-9,
                  Polarization = 0,
                  E0           = 1)

Detector = Photodiode(NA                = 0.8,
                     Sampling          = 1001,
                     GammaOffset       = 0,
                     PhiOffset         = 0)


Detector.Plot()
