
"""
_________________________________________________________
Optimization of RI dependence minimizing STD of detector response.
_________________________________________________________
"""



import matplotlib.pyplot as plt
import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.functions.Optimization import OptimizeRI

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
            PhiOffset   = 15
            )

LP01.magnificate(magnification=2.)

DiameterList = np.linspace(100,1000,10).round(4) * 1e-9

RIList = np.linspace(1.3, 2.0, 4).round(4)

SourceKwargs = {'wavelength': 400e-9,
                'npts': 101,
                'Meshes': LP01.Meshes}

DataFrame = OptimizeRI(RIList,
                       DiameterList,
                       Detector = LP01,
                       **SourceKwargs)

DataFrame.xs('Parallel').unstack(1).plot(y='Values',
                                         figsize=(8,3),
                                         grid=True,
                                         title='Parallel coupling')

DataFrame.xs('Parallel').unstack(1).plot(y='STD',
                                         figsize=(8,3),
                                         grid=True,
                                         title='Parallel coupling STD')

DataFrame.xs('Perpendicular').unstack(1).plot(y='Values',
                                              figsize=(8,3),
                                              grid=True,
                                              title='Parallel coupling')

DataFrame.xs('Perpendicular').unstack(1).plot(y='STD',
                                              figsize=(8,3),
                                              grid=True,
                                              title='Perpendicular coupling')

plt.show()















# -













# -
