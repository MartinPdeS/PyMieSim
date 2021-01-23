
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Optimizer import Simulator
from PyMieCoupling.classes.Sets import ExperimentalSet

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)



Photodiode0 = Photodiode(NA                = 0.5,
                         Source            = LightSource,
                         Sampling          = 801,
                         GammaOffset       = 0,
                         PhiOffset         = 0,
                         CouplingMode      = 'Mean')


Set = ExperimentalSet(DiameterList  = np.linspace(100,1000,100).round(4) * 1e-9,
                      RIList        = np.linspace(1.3, 1.5, 4).round(4),
                      Detectors     = Photodiode0,
                      Source        = LightSource,)

def EvalFunc(x):

    Set.Detectors['Detector 0'].NA = x[0]

    Array = Set.Coupling()

    return Array.Cost('Max') # can be: RI_STD  -  RI_RSD  -  Monotonic  -  Mean  -  Max  -  Min


Minimizer = Simulator(EvalFunc, ParameterName= ['NA'])

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [0.2],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-5,
                  options  = {'maxiter': 50, 'rhobeg':0.1})

print(Result)

DF = Set.GetCouplingFrame()

DF.Plot('Coupling') # can be Couplimg  -  STD

plt.show()







# -
