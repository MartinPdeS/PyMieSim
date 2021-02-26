from PyMieSim.Source import PlaneWave
from PyMieSim.Detector import LPmode
from PyMieSim.Scatterer import Sphere

LightSource = PlaneWave(Wavelength = 450e-9,
                        Polarization = 0)

Detector = LPmode(Mode         = (1, 1,'h'),
                  Sampling     = 201,
                  NA           = 0.2,
                  GammaOffset  = 0,
                  PhiOffset    = 0,
                  CouplingMode = 'Centered')


Scat = Sphere(Diameter    = 300e-9,
              Source      = LightSource,
              Index       = 1.4)

Coupling = Detector.Coupling(Scatterer = Scat)

print(Coupling)
