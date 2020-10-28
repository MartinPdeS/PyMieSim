
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""



import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import sys
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.functions.Optimization import CouplingStat



Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)


DiameterList = np.linspace(100,1000,10).round(4) * 1e-9

RIList = np.linspace(1.3, 2.0, 10).round(4)

def EvalFunc(x):

    LP01 = LPmode(Fiber       = Fiber,
                  Mode        = (0, 1),
                  Wavelength  = 400e-9,
                  Npts        = 101,
                  ThetaOffset = 0,
                  PhiOffset   = x)

    DataFrame = CouplingStat(RIList       = RIList,
                             DiameterList = DiameterList,
                             Detector     = LP01,
                             Wavelength   = 400e-9)

    print('\n-> PhiOffset: {0}\n-> Max coupling: {1}\n'.format(x, DataFrame.ParaMax), flush=True)

    return DataFrame.ParaMax


Result = optimize.minimize_scalar(EvalFunc, options={'maxiter':5})

XMin = Result.x












# -













# -
