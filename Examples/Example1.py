import matplotlib.pyplot as plt
import numpy as np
from progress.bar import Bar

from miecoupling.src.classes.Fiber import fiber
from miecoupling.src.classes.Modes import mode
from miecoupling.src.classes.Scattering import Scatterer
from miecoupling.src.functions.couplings import PointFieldCoupling, MeanFieldCoupling

npts=101

Fiber = fiber(4.2e-6,
              1.4456,
              20.5e-6,
              1.4444)

LP11 = mode(fiber=Fiber,
            LPmode=(2, 1),
            wavelength=400e-9,
            npts=npts,
            )

LP11.magnificate(magnification=2.)


LP01 = mode(fiber=Fiber,
            LPmode=(1, 1),
            wavelength=400e-9,
            npts=npts,
            )

LP01.magnificate(magnification=2.)





DiameterList = np.linspace(100,9000,50) * 1e-9



ThetaOffset, PhiOffset = 0,+10

CouplingLP01 = []
CouplingLP11 = []
with Bar('Processing...', max = len(DiameterList)) as bar:
    for Diameter in DiameterList:

        Scat = Scatterer(diameter=Diameter,
                         wavelength=400e-9,
                         index=1.4,
                         npts=400)

        Scat.GenField(PolarizationAngle=0)

        Scat.SampleField(ThetaBound=[-10+ThetaOffset,10+ThetaOffset],
                         PhiBound=[-10+PhiOffset, 10+PhiOffset],
                         npts=npts)

        CouplingLP01.append( PointFieldCoupling(LP01.Fourier, Scat.Sample.Perpendicular) )
        CouplingLP11.append( PointFieldCoupling(LP11.Fourier, Scat.Sample.Perpendicular) )
        bar.next()


fig = plt.figure(figsize=(15,5))
ax0 = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

ax0.plot(DiameterList*1e6, CouplingLP01, 'C0', label=r'LP$_{01}$')
ax0.plot(DiameterList*1e6, CouplingLP11, 'C1', label=r'LP$_{11+}$')
ax0.set_title('Mode coupling vs. Scatterer diameter')
ax0.legend()
ax0.set_xlabel(r'Scatterer diameter [$\mu m$]')
ax0.set_ylabel('Modal Coupling')
ax0.set_yscale('log')
ax0.grid()

ax1.plot(DiameterList*1e6, np.array(CouplingLP11)/np.array(CouplingLP01))
ax1.set_title('Coupling ratio')
ax1.set_yscale('log')
ax1.grid()



plt.show()










# -
