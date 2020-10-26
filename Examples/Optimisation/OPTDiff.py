
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
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.functions.Optimization import CouplingStat

npts=51

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)


LP01 = mode(fiber       = Fiber,
            LPmode      = (0, 1),
            wavelength  = 400e-9,
            npts        = npts,
            ThetaOffset = 0,
            PhiOffset   = 35
            )

LP01.magnificate(magnification=2.)

DiameterList = np.linspace(100,1000,10).round(4) * 1e-9

RIList = np.linspace(1.3, 2.0, 10).round(4)

SourceKwargs = {'wavelength': 400e-9,
                'npts': 101,
                'Meshes': LP01.Meshes}


def EvalFunc(x):

    LP01 = mode(fiber       = Fiber,
                LPmode      = (0, 1),
                wavelength  = 400e-9,
                npts        = npts,
                ThetaOffset = 0,
                PhiOffset   = x
                )

    DataFrame = CouplingStat(RIList,
                             DiameterList,
                             Detector = LP01,
                             **SourceKwargs)

    print('\n-> PhiOffset: {0}\n-> Max coupling: {1}\n'.format(x, DataFrame.ParaMax), flush=True)

    return DataFrame.ParaMax


Result = optimize.minimize_scalar(EvalFunc, options={'maxiter':5})

XMin = Result.x












# -













# -
