import numpy as np
from PyMieSim            import Material
from PyMieSim.Scatterer  import Sphere
from PyMieSim.Source     import PlaneWave
from PyMieSim.Detector   import Photodiode
from PyMieSim.Experiment import ScatSet, SourceSet, Setup, DetectorSet

scatKwargs   = { 'Diameter'    : np.linspace(400e-9, 2000e-9, 200),
                 'Material'    : Material('BK7'),
                 'nMedium'     : [1] }

sourceKwargs = { 'Wavelength'   : 1e-6,
                 'Polarization' : [0]}

Detector0 = Photodiode(NA                = 2.0,
                       Sampling          = 300,
                       GammaOffset       = 0,
                       PhiOffset         = 0,
                       CouplingMode      = 'Centered')

detecSet   = DetectorSet([Detector0])

scatSet    = ScatSet(Scatterer = Sphere,  kwargs = scatKwargs )

sourceSet  = SourceSet(Source = PlaneWave, kwargs = sourceKwargs )

Experiment = Setup(ScattererSet = scatSet,
                   SourceSet    = sourceSet,
                   DetectorSet  = detecSet)

Coupling = Experiment.Coupling(AsType='pymiesim')

Coupling.Plot(x='Diameter')
