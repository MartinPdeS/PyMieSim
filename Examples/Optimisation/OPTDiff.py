
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""

import numpy as np
from scipy import optimize
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import LPmode
from PyMieCoupling.functions.Optimization import CouplingStat
from PyMieCoupling.classes.Misc import Source

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP01 = LPmode(Fiber       = Fiber,
              Mode        = (0, 1),
              Source      = LightSource,
              Npts        = 201,
              ThetaOffset = 0,
              PhiOffset   = 0,
              GPU         = False)



def EvalFunc(x):

    LP01.PhiOffset = x

    DataFrame = CouplingStat(RIList       = np.linspace(1.3, 2.0, 3).round(4),
                             DiameterList = np.linspace(100,1000,3).round(4) * 1e-9,
                             Detector     = LP01,
                             Source       = LightSource,
                             GPU          = False)

    print('\n-> PhiOffset:    {0}\n-> Max coupling: {1}\n'.format(x, DataFrame.ParaMax), flush=True)

    return DataFrame.ParaMax


Result = optimize.minimize_scalar(EvalFunc, options={'maxiter':5})

XMin = Result.x












# -













# -
