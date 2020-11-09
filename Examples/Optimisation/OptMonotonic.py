
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

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP01 = LPmode(Fiber       = Fiber,
              Mode        = (1, 1),
              Source      = LightSource,
              Npts        = 51,
              ThetaOffset = 0,
              PhiOffset   = 0,
              cuda        = False)



RIList = [1.5]
DiameterList = np.linspace(100,1000,100).round(4) * 1e-9



def EvalFunc(X):

    LP01.PhiOffset = X[0]

    LP01.ThetaOffset = X[1]

    Array = LoopRIDiameter(RIList       = RIList,
                           DiameterList = DiameterList,
                           Detector     = LP01,
                           Source       = LightSource,
                           cuda         = False)


    return Array.Monotonic()


Minimizer = Simulator(EvalFunc)

Result = minimize(fun      = Minimizer.simulate,
                  x0       = [30, 30],
                  method   = 'COBYLA',
                  callback = Minimizer.callback,
                  tol      = 1e-10,
                  options  = {'maxiter': 20, 'rhobeg':20})


print(Result)

LP01.PhiOffset = Result.x[0]

LP01.ThetaOffset = Result.x[1]

DF = Frame(RIList        = RIList,
           DiameterList  = DiameterList,
           Detector      = LP01,
           Source        = LightSource,
           cuda          = False)

DF.plot(y='Coupling')


plt.show()





# -
