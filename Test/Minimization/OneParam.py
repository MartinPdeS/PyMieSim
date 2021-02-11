
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
from scipy.optimize import minimize
from PyMieSim.classes.Detector import Photodiode, LPmode
from PyMieSim.Physics import Source
from PyMieSim.classes.Optimizer import Simulator
from PyMieSim.classes.Sets import ExperimentalSet, ScattererSet

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)



Detector0 = Photodiode(NA                = 0.2,
                       Sampling          = 150,
                       GammaOffset       = 0,
                       PhiOffset         = 0,
                       CouplingMode      = 'Centered')

Detector1 = LPmode(NA                = 0.2,
                   Sampling          = 150,
                   Mode              = (0,1),
                   GammaOffset       = 0,
                   PhiOffset         = 0,
                   CouplingMode      = 'Centered')


ScatSet = ScattererSet(DiameterList  = np.linspace(100e-9, 3500e-9, 100),
                       RIList        = np.linspace(1.5, 1.5, 1).round(1),
                       Source        = LightSource)

Set = ExperimentalSet(ScattererSet  = ScatSet,
                      Detectors     = (Detector0, Detector1))


def EvalFunc(x):

    Set.Detectors[1].NA = x[0]

    return Set.Coupling.Cost('Max') # can be: RI_STD  -  RI_RSD  -  Monotonic  -  Mean  -  Max  -  Min


Minimizer = Simulator(EvalFunc, ParameterName= ['NA'])

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [0.2],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-5,
                  options  = {'maxiter': 10, 'rhobeg':0.1})

print(Result)

Set.DataFrame.Plot('Coupling') # can be Couplimg  -  STD








# -
