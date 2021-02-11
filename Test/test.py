from PyMieSim.Physics import Source
from PyMieSim.classes.Detector import LPmode
from PyMieSim.classes.Scattering import Scatterer

LightSource = Source(Wavelength = 450e-9,
                  Polarization = 0,
                  Power = 1,
                  Radius = 1)

Detector = LPmode(Mode         = (1, 1,'h'),
                  Sampling     = 201,
                  NA           = 0.2,
                  GammaOffset  = 0,
                  PhiOffset    = 0,
                  CouplingMode = 'Centered')


Scat = Scatterer(Diameter    = 400e-9,
                 Source      = LightSource,
                 Index       = 1.4)

Coupling = Detector.Coupling(Scatterer = Scat)
