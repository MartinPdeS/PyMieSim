from numpy           import linspace
from PyMieSim.Source import PlaneWave
from PyMieSim.Sets   import ScattererSet

DiameterList  = linspace(100e-9, 10000e-9, 400)
IndexList     = linspace(1.5, 1.8, 3).round(1)

Source = PlaneWave(Wavelength   = 450e-9,
                  Polarization = 0,
                  E0           = 1)


ScatSet = ScattererSet(DiameterList = DiameterList,
                       IndexList    = IndexList,
                       Source       = Source)


Qsca = ScatSet.Qsca()

fig = Qsca.Plot()
