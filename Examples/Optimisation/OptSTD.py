
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.functions.Optimization import ComputeSTD
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.DataFrame import CouplingStat

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP01 = LPmode(Fiber       = Fiber,
              Mode        = (0, 1),
              Source      = LightSource,
              Npts        = 51,
              ThetaOffset = 0,
              PhiOffset   = 0,
              cuda        = False)



RIList = np.linspace(1.3, 2.0, 6).round(4)
DiameterList = np.linspace(100,1000,50).round(4) * 1e-9

def EvalFunc(x):

    LP01.PhiOffset = x

    STD = ComputeSTD(RIList       = RIList,
                     DiameterList = DiameterList,
                     Detector     = LP01,
                     Source       = LightSource,
                     cuda         = False)

    print('\n-> PhiOffset:    {0}\n-> Max coupling: {1}\n'.format(x, np.sum(STD)), flush=True)

    return np.sum(STD)


Result = optimize.minimize_scalar(EvalFunc, options={'maxiter':15})


LP01.PhiOffset = Result.x

DataFrame = CouplingStat(RIList        = RIList,
                         DiameterList  = DiameterList,
                         Detector      = LP01,
                         Source        = LightSource,
                         cuda          = False)

DataFrame.plot(y='Coupling')

DataFrame.plot(y='STD')

plt.show()









# -













# -
