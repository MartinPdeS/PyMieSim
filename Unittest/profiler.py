
"""
_________________________________________________________
Profiler code for optimization.
_________________________________________________________
"""



import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling

npts=201

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = mode(fiber      = Fiber,
            LPmode     = (2, 1),
            wavelength = 400e-9,
            npts       = npts,
            )

LP01 = mode(fiber       = Fiber,
            LPmode      = (1, 1),
            wavelength  = 400e-9,
            npts        = npts,
            )


LP01.magnificate(magnification=2.)

LP11.magnificate(magnification=2.)

DiameterList = np.linspace(100,1000,10) * 1e-9

CouplingLP01, CouplingLP11 = [], []
def main():
    for Diameter in tqdm(DiameterList, total = len(DiameterList), desc ="Progress"):

        Scat = Scatterer(diameter    = Diameter,
                         wavelength  = 400e-9,
                         index       = 1.4,
                         npts        = npts,
                         Meshes      = LP01.Meshes,
                         CacheTrunk  = None)

        CouplingLP01.append( PointFieldCoupling(Detector = LP01,
                                                Source   = Scat.Field.Parallel,
                                                Mesh     = Scat.Meshes) )


from pycallgraph import PyCallGraph
from pycallgraph.output import GraphvizOutput
with PyCallGraph(output=GraphvizOutput()):
    main()













# -
