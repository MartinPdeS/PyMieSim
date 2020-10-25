
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

LP01Perp, LP01Para, LP11Perp, LP11Para = [], [], [], []

for Diameter in tqdm(DiameterList, total = len(DiameterList), desc ="Progress"):

    Scat = Scatterer(diameter    = Diameter,
                     wavelength  = 400e-9,
                     index       = 1.4,
                     Meshes      = LP01.Meshes
                     )
    CLP01 = PointFieldCoupling(Detector = LP01, Source   = Scat)
    CLP11 = PointFieldCoupling(Detector = LP11, Source   = Scat)


    LP01Perp.append( CLP01[0] )
    LP01Para.append( CLP01[1] )
    LP11Perp.append( CLP11[0] )
    LP11Para.append( CLP11[1] )





fig, (ax0, ax1) = plt.subplots(2,1, figsize=(15,8), sharex = True, subplot_kw={'title':'Mode coupling vs. Scatterer diameter'})

ax0.plot(DiameterList*1e6, LP01Perp, 'C0', label=r'LP$_{01}$ Perp.')
ax0.plot(DiameterList*1e6, LP01Para, 'C1', label=r'LP$_{01}$ Para.')

ax0.plot(DiameterList*1e6, LP11Perp, 'C2', label=r'LP$_{11+}$ Perp.')
ax0.plot(DiameterList*1e6, LP11Para, 'C3', label=r'LP$_{11+}$ Para.')

ax0.legend()

ax1.set_xlabel(r'Scatterer diameter [$\mu m$]')

ax0.set_ylabel('Modal Coupling')

ax0.set_yscale('log')

ax0.grid()

ax1.plot(DiameterList*1e6, np.array(LP11Perp)/np.array(LP01Perp), 'C0', label='Perp')
ax1.plot(DiameterList*1e6, np.array(LP11Para)/np.array(LP01Para), 'C1', label='Para')

ax1.set_ylabel(r'Modal coupling ratio $\Gamma_{11+}$')

ax1.legend()

ax1.grid()

plt.show()









# -
