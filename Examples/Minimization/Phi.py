
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.classes.Fields import Source
from PyMieCoupling.classes.Optimizer import Simulator1 as Simulator
from PyMieCoupling.classes.Scattering import ScattererSet

LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)


LP01 = LPmode(Mode          = (0, 1),
               Source        = LightSource,
               Npts          = 81,
               ThetaOffset   = 0,
               PhiOffset     = 0,
               Name          = 'LP01',
               NA            = 0.9)


LP11 = LPmode(Mode          = (1, 1),
              Source        = LightSource,
              Npts          = 81,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Name          = 'LP11',
              NA            = 0.54)


Photodiode0 = Photodiode(NA                = 0.54,
                         Source            = LightSource,
                         Npts              = 81,
                         ThetaOffset       = 0,
                         PhiOffset         = 0)


Set = ScattererSet(DiameterList  = np.linspace(100,1000,100).round(4) * 1e-9,
                   RIList        = np.linspace(1.3, 1.5, 6).round(4),
                   Detector      = LP01,
                   Source        = LightSource,
                   Mode          = 'Mean'
                   )

def EvalFunc(x):

    Set.Detector.PhiOffset = x

    Array = Set.GetCoupling(Filter = 'None')

    return Array.Cost('Max') # can be: RI_STD  -  RI_RSD  -  Monotonic  -  Mean  -  Max  -  Min


Minimizer = Simulator(EvalFunc)

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [0],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-5,
                  options  = {'maxiter': 50, 'rhobeg':20})
print(Result)

DF = Set.GetFrame()

DF.Plot('Coupling') # can be Couplimg  -  STD

plt.show()







# -
