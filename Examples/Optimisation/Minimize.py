
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
              Npts        = 81,
              ThetaOffset = 0,
              PhiOffset   = 0,
              cuda        = False,
              Magnification = 0.1)



RIList = np.linspace(1.3, 1.6, 1).round(4)
DiameterList = np.linspace(10,100,100).round(4) * 1e-9

def EvalFunc(x):

    LP01.PhiOffset = x[0]

    LP01.ThetaOffset = x[1]

    Array = LoopRIDiameter(RIList       = RIList,
                           DiameterList = DiameterList,
                           Detector     = LP01,
                           Source       = LightSource,
                           cuda         = False,
                           QuietMode    = True,
                           Polarization = 'Perpendicular')


    return Array.Cost('Mean') # can be: RI  -  RI/Mean  -  Monotonic  -  Mean  -  max


Minimizer = Simulator(EvalFunc)

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [10, 10],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-3,
                  options  = {'maxiter': 50, 'rhobeg':50})
print(Result)

#LP01.PhiOffset = Result.x[0]

LP01.ThetaOffset = Result.x[1]

print(LP01.Meshes.Phi.Boundary.Degree)
DF = Frame(RIList        = RIList,
           DiameterList  = DiameterList,
           Detector      = LP01,
           Source        = LightSource,
           cuda          = False)

DF.plot(y='Coupling')

DF.plot(y='STD')

plt.show()









# -
