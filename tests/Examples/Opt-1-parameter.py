import numpy as np
from PyMieSim.Detector import Photodiode, LPmode
from PyMieSim.Source import PlaneWave
from PyMieSim.Optimization import Optimizer
from PyMieSim.Experiment import ScatSet, SourceSet, Setup

DiameterList   = np.linspace(100e-9, 1000e-9, 200)

Detector0 = Photodiode(NA                 = 0.1,
                     Sampling          = 300,
                     GammaOffset       = 20,
                     PhiOffset         = 0,
                     CouplingMode      = 'Centered')

scat = ScatSet(DiameterList   = DiameterList,
             IndexList      = [1.5],
             nMedium        = 1,
             ScattererType  = 'Sphere')

source = SourceSet(WavelengthList  = 400e-9,
                 PolarizationList  = [0],
                 SourceType        = 'PlaneWave')


Experiment = Setup(ScattererSet = scat,
                 SourceSet    = source,
                 DetectorSet  = [Detector0])


# Metric can be "max" - "min" - "mean"
#"std+RI" - "std+Diameter" - "std+Polarization" - "std+Wavelength" - "std+Detector"
#"monotonic+RI" - "monotonic+Diameter" - "monotonic+Polarization" - "monotonic+Wavelength" - "monotonic+Detector"

Opt    = Optimizer(Setup           = Experiment,
                 Metric          = 'mean',
                 Parameter       = ['PhiOffset'],
                 Optimum         = 'Maximum',
                 MinVal          = [1e-5],
                 MaxVal          = [180],
                 WhichDetector   = 0,
                 X0              = [0.6],
                 MaxIter         = 350,
                 Tol             = 1e-4,
                 FirstStride     = 30)

print(Opt.Result)

df = Experiment.Coupling(AsType='dataframe')

df.Plot(y='Coupling', x='Diameter') # can be "Couplimg"  or  "STD"
