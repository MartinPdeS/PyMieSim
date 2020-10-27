
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import matplotlib.pyplot as plt
import numpy as np
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Detector import mode
from PyMieCoupling.classes.Scattering import Scatterer
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
            Magnification = 2.)

DiameterList = np.linspace(100,1000,10).round(3) * 1e-9


def main():

    DataFrame = CouplingStat(RIList       = [1.4],
                             DiameterList = DiameterList,
                             Detector     = LP01,
                             wavelength   = 400e-9,
                             npts         = 101,
                             debugMode    = False)


from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
with PyCallGraph(output=GraphvizOutput()):
    main()













# -
