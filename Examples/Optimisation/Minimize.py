
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.functions.Optimization import LoopRIDiameter
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.DataFrame import Frame
from PyMieCoupling.classes.Optimizer import Simulator

LightSource = Source(Wavelength   = 1000e-9,
                     Polarization = 0)

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP01 = LPmode(Fiber       = Fiber,
              Mode        = (0, 1),
              Source      = LightSource,
              Npts        = 101,
              ThetaOffset = 0,
              PhiOffset   = 0,
              Magnification = 0.1)

RIList = np.linspace(1.3, 1.6, 4).round(4)

DiameterList = np.linspace(10,400,50).round(4) * 1e-9

def EvalFunc(x):

    LP01.PhiOffset = x[0]

    LP01.ThetaOffset = x[1]

    Array = LoopRIDiameter(RIList       = RIList,
                           DiameterList = DiameterList,
                           Detector     = LP01,
                           Source       = LightSource,
                           QuietMode    = True,
                           Polarization = 'Perpendicular')


    return Array.Cost('RI') # can be: RI  -  RI/Mean  -  Monotonic  -  Mean  -  Max  -  Min


Minimizer = Simulator(EvalFunc)

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [10, 10],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-20,
                  options  = {'maxiter': 150, 'rhobeg':50})
print(Result)

LP01.PhiOffset = Result.x[0]

LP01.ThetaOffset = Result.x[1]



DF = Frame(RIList        = RIList,
           DiameterList  = DiameterList,
           Detector      = LP01,
           Source        = LightSource)

DF.plot(y='Coupling')

DF.plot(y='STD')

plt.show()









# -
