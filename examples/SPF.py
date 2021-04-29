from PyMieSim.Scatterer import Sphere
from PyMieSim.Source    import PlaneWave

Source = PlaneWave(Wavelength   = 450e-9,
                  Polarization = 0,
                  E0           = 1)

Scat = Sphere(Diameter    = 800e-9,
            Source      = Source,
            Index       = 1.4)

SPF = Scat.SPF(Num=100)

SPF.Plot()
