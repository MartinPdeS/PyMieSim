
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Optimizer import Simulator
from PyMieCoupling.classes.Scattering import ScattererSet

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)


LP01 = LPmode(Name        = 'LP11',
              Mode        = (0, 1),
              Source      = LightSource,
              Samples     = 501,
              NA          = 0.5,
              GammaOffset = 0,
              PhiOffset   = 0
              )


LP11 = LPmode(Name        = 'LP11',
              Mode        = (1, 1),
              Source      = LightSource,
              Samples     = 501,
              NA          = 0.5,
              GammaOffset = 0,
              PhiOffset   = 0
              )


Photodiode0 = Photodiode(Name              = 'Detector',
                      NA                = 0.5,
                      Source            = LightSource,
                      Samples           = 501,
                      GammaOffset       = 0,
                      PhiOffset         = 0)


Set = ScattererSet(DiameterList  = np.linspace(100,1000,100).round(4) * 1e-9,
                   RIList        = np.linspace(1.3, 1.5, 6).round(4),
                   Detector      = Photodiode0,
                   Source        = LightSource,
                   Mode          = 'Mean'
                   )

def EvalFunc(x):

    Set.Detector.Phi = x

    Array = Set.GetCouplingArray()

    return Array.Cost('Max') # can be: RI_STD  -  RI_RSD  -  Monotonic  -  Mean  -  Max  -  Min


Minimizer = Simulator(EvalFunc, ParameterName= ['Phi'])

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [0],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-5,
                  options  = {'maxiter': 50, 'rhobeg':20})
print(Result)

DF = Set.GetCouplingFrame()

DF.Plot('Coupling') # can be Couplimg  -  STD

plt.show()







# -
