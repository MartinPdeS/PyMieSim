
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""



import matplotlib.pyplot as plt
import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.functions.Optimization import CouplingStat

npts=201

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

DiameterList = np.linspace(100,1000,20).round(2) * 1e-9

RIList = np.linspace(1.3, 2.0, 10).round(2)

SourceKwargs = {'wavelength': 400e-9,
                'npts': 101}

DataFrame = CouplingStat(RIList,
                         DiameterList,
                         Detector = LP01,
                         **SourceKwargs)

DataFrame.plot(y='Coupling')

DataFrame.plot(y='STD')

plt.show()















# -













# -
