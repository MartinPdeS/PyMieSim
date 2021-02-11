
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""

import matplotlib.pyplot as plt
import numpy as np
from PyMieSim.classes.Fiber import fiber
from PyMieSim.classes.Detector import LPmode
from PyMieSim.classes.Scattering import Scatterer
from PyMieSim.classes.Misc import Source
from PyMieSim.classes.DataFrame import Frame
from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
#from PyMieSim.functions.MieComputing import MieS1S2
from PyMieSim.S1S2 import MieS1S2

npts = 51
cuda = False

LightSource = Source(Wavelength   = 400e-9,
                     Polarization = 0)


Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)


Detector0 = LPmode(Fiber         = Fiber,
                   Mode          = (0, 1),
                   Source        = LightSource,
                   Npts          = npts,
                   ThetaOffset   = 0,
                   PhiOffset     = 35,
                   Magnification = 2.,
                   cuda          = cuda)

DiameterList = np.linspace(100,1000,10).round(3) * 1e-9


def main():
    for i in range(100):
        MieS1S2(1.4, 10, 0.5)




with PyCallGraph(output=GraphvizOutput()):
    main()













# -
