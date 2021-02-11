
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieSim.classes.Detector import LPmode, Photodiode
from PyMieSim.classes.Fields import Source
from PyMieSim.classes.Optimizer import Simulator
from PyMieSim.classes.Scattering import ScattererSet

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
              Npts          = 41,
              ThetaOffset   = 0,
              PhiOffset     = 0,
              Name          = 'LP11',
              NA            = 0.2,
              Orientation   = 'v')


Photodiode0 = Photodiode(NA                = 0.54,
                         Source            = LightSource,
                         Npts              = 41,
                         ThetaOffset       = 0,
                         PhiOffset         = 0)


Set = ScattererSet(DiameterList  = np.linspace(100,1000,100).round(4) * 1e-9,
                   RIList        = np.linspace(1.3, 1.5, 6).round(4),
                   Detector      = LP11,
                   Source        = LightSource,
                   Mode          = 'Centered'
                   )

def EvalFunc(x):

    Set.Detector.PhiOffset = x[0]

    Set.Detector.NA = x[1]

    Array = Set.GetCouplingArray(Filter = 'None')

    return Array.Cost('Max') # can be: RI_STD  -  RI_RSD  -  Monotonic  -  Mean  -  Max  -  Min


Minimizer = Simulator(EvalFunc, ParameterName= ['Phi', 'NA'])

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [0, 0.2],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-5,
                  options  = {'maxiter': 50, 'rhobeg':0.5})
print(Result)

DF = Set.GetCouplingFrame()

DF.Plot('Coupling') # can be Couplimg  -  STD

plt.show()







# -
