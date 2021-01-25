
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""



from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Scattering import Scatterer

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

LP11 = LPmode(Mode         = (1, 1),
              Sampling     = 501,
              NA           = 0.2,
              GammaOffset  = 0,
              PhiOffset    = 90,
              CouplingMode = 'Centered')

Scat = Scatterer(Diameter    = 50e-9,
                 Source      = LightSource,
                 Index       = 1.4)


print(LP11.Coupling(Scat))


# -
