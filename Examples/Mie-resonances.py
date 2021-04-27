from PyMieSim.Scatterer import Sphere
from PyMieSim.Detector import Photodiode
from PyMieSim.Experiment import ScatSet, SourceSet, Setup
import numpy as np

WavelengthList   = np.linspace(400e-9, 1000e-9, 200)

Detector0 = Photodiode(NA               = 0.1,
                     Sampling          = 300,
                     GammaOffset       = 20,
                     PhiOffset         = 30,
                     CouplingMode      = 'Centered')

scat = ScatSet(DiameterList   = 200e-9,
              IndexList         = [1.5],
              nMedium        = 1,
              ScattererType  = 'Sphere')

source = SourceSet(WavelengthList   = WavelengthList,
                 PolarizationList  = [0],
                 SourceType        = 'PlaneWave')


Experiment = Setup(ScattererSet = scat,
                  SourceSet    = source,
                  DetectorSet  = [Detector0])

Qsca = Experiment.Qsca(AsType='dataframe')

Qsca.Plot(y='Qsca', x='Wavelength')
