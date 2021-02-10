from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.Physics import Source

LightSource = Source(Wavelength = 450e-9,
                  Polarization = 0,
                  Power = 1,
                  Radius = 1)

Scat = Scatterer(Diameter    = 300e-9,
                 Source      = LightSource,
                 Index       = 1.4)


Fields = Scat.Field(Num=100)

Fields.Plot()
