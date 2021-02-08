from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.Scattering import Scatterer

LightSource = Source(Wavelength = 450e-9,
                  Polarization = 0,
                  Power = 1,
                  Radius = 1)

Detector0 = Photodiode(NA               = 0.2,
                      Sampling          = 201,
                      GammaOffset       = 0,
                      PhiOffset         = 0,
                      CouplingMode = 'Centered')

Detector = LPmode(Mode         = (0, 1,'v'),
                  Sampling     = 501,
                  NA           = 0.2,
                  GammaOffset  = 0,
                  PhiOffset    = 0,
                  CouplingMode = 'Centered')


Scat = Scatterer(Diameter    = 1e-9,
                 Source      = LightSource,
                 Index       = 1.4)



Detector.Coupling(Scat)
#Scat.Field().Plot()
#Footprint = Detector.Footprint(Scatterer = Scat)


#Footprint.Plot()
#Footprint
