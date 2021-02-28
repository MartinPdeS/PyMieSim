from PyMieSim.Source import PlaneWave
from PyMieSim.Detector import Photodiode

LightSource = PlaneWave(Wavelength = 450e-9,
                        Polarization = 0)

Detector = Photodiode(NA                = 0.8,
                      Sampling          = 1001,
                      GammaOffset       = 0,
                      PhiOffset         = 0)


Detector.Plot()
