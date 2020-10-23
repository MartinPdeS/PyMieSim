
"""
_________________________________________________________
Scattering Parallel Field coupling with en-face LP01 and LP11 Mode
For different scatterer diameters.
_________________________________________________________
"""



import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from PyMieCoupling.classes.Fiber import fiber
from PyMieCoupling.classes.Modes import mode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling
from PyMieCoupling.functions.Optimization import OptimizeRI

npts=201

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

LP11 = mode(fiber       = Fiber,
            LPmode      = (1, 1),
            wavelength  = 400e-9,
            npts        = npts,
            ThetaOffset = 0,
            PhiOffset   = 0
            )

LP01 = mode(fiber       = Fiber,
            LPmode      = (0, 1),
            wavelength  = 400e-9,
            npts        = npts,
            ThetaOffset = 0,
            PhiOffset   = 0
            )



LP01.magnificate(magnification=2.)

LP11.magnificate(magnification=2.)

DiameterList = np.linspace(100,1000,2) * 1e-9

CouplingLP01, CouplingLP11 = [], []

for Diameter in tqdm(DiameterList, total = len(DiameterList), desc ="Progress"):

    Scat = Scatterer(diameter    = Diameter,
                     wavelength  = 400e-9,
                     index       = 1.4,
                     Meshes      = LP01.Meshes
                     )

    CouplingLP01.append( PointFieldCoupling(Detector = LP01,
                                            Source   = Scat,
                                            Field = 'Parallel') )

    CouplingLP11.append( PointFieldCoupling(Detector = LP11,
                                            Source   = Scat,
                                            Field = 'Parallel') )




fig, (ax0, ax1) = plt.subplots(2,1, figsize=(15,8), sharex = True, subplot_kw={'title':'Mode coupling vs. Scatterer diameter'})

ax0.plot(DiameterList*1e6, CouplingLP01, 'C0', label=r'LP$_{01}$')

ax0.plot(DiameterList*1e6, CouplingLP11, 'C1', label=r'LP$_{11+}$')

ax0.legend()

ax1.set_xlabel(r'Scatterer diameter [$\mu m$]')

ax0.set_ylabel('Modal Coupling')

ax0.set_yscale('log')

ax0.grid()

ax1.plot(DiameterList*1e6, np.array(CouplingLP11)/np.array(CouplingLP01))

ax1.set_ylabel(r'Modal coupling ratio $\Gamma_{11+}$')

ax1.grid()

plt.show()









# -
