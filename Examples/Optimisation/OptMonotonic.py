
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
from PyMieCoupling.functions.Optimization import ComputeMonotonic
from PyMieCoupling.classes.Misc import Source
from PyMieCoupling.classes.DataFrame import CouplingStat
from PyMieCoupling.classes.Optimizer import Simulator

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



RIList = [1.5]
DiameterList = np.linspace(100,1000,100).round(4) * 1e-9



def EvalFunc(Angle):

    LP01.PhiOffset = Angle[0]

    LP01.ThetaOffset = Angle[1]

    val = ComputeMonotonic(RIList       = RIList,
                           DiameterList = DiameterList,
                           Detector     = LP01,
                           Source       = LightSource,
                           cuda         = False,
                           QuietMode    = True)

    return val


Minimizer = Simulator(EvalFunc)

Result = minimize(Minimizer.simulate,
                  [30, 30],
                  method='COBYLA',
                  callback=Minimizer.callback,
                  tol    = 1e-2,
                  options={"disp": True, 'rhobeg':10})


LP01.PhiOffset = Result.x[0]

LP01.ThetaOffset = Result.x[1]

DataFrame = CouplingStat(RIList        = RIList,
                         DiameterList  = DiameterList,
                         Detector      = LP01,
                         Source        = LightSource,
                         cuda          = False)

DataFrame.plot(y='Coupling')

DataFrame.plot(y='STD')

plt.show()





# -
