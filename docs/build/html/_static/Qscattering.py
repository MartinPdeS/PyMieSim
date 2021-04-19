import numpy as np
from PyMieSim.Source import PlaneWave
from PyMieSim.Sets import ScattererSet


Source = PlaneWave(Wavelength   = 450e-9,
                  Polarization = 0,
                  E0           = 1)


ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 10000e-9, 400),
                      IndexList        = np.linspace(1.5, 1.8, 3).round(1),
                      Source        = Source)


Qsca = ScatSet.Qsca()

fig = Qsca.Plot()
