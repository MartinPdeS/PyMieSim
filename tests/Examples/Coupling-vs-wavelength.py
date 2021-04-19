from PyMieSim.Scatterer import Sphere
from PyMieSim.Detector import Photodiode
from PyMieSim.Experiment import ScatSet, SourceSet, Setup
import numpy as np

WavelengthList = np.linspace(400e-9, 1000e-9, 100)

Detector0 = Photodiode(NA                = 2.0,
                     Sampling          = 300,
                     GammaOffset       = 0,
                     PhiOffset         = 0,
                     CouplingMode      = 'Centered')

scat = ScatSet(DiameterList  = [200e-9],
             IndexList        = [4],
             nMedium       = 1,
             ScattererType = 'Sphere')

source = SourceSet(WavelengthList   = WavelengthList,
                 PolarizationList = [0],
                 SourceType       = 'PlaneWave')


Experiment = Setup(ScattererSet = scat,
                 SourceSet    = source,
                 DetectorSet  = [Detector0])

DF = Experiment.Coupling(AsType='dataframe')

DF.Plot(y='Coupling', x='Wavelength')
