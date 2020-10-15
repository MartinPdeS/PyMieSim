import matplotlib.pyplot as plt
import numpy as np
from progress.bar import Bar

from src.classes.Fiber import fiber
from src.classes.Modes import mode
from src.classes.Scattering import Scatterer
from src.functions.couplings import PointFieldCoupling, MeanFieldCoupling

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





DiameterList = np.linspace(100,900,150) * 1e-9

print(LP01.AngleBound)
CouplingLP01 = []
CouplingLP11 = []
with Bar('Processing...', max = len(DiameterList)) as bar:
    for Diameter in DiameterList:

        Scat = Scatterer(diameter=Diameter,
                         wavelength=400e-9,
                         index=1.4,
                         npts=400)

        Scat.GenField(PolarizationAngle=10)

        Scat.SampleField(ThetaBound=[180,190],
                         PhiBound=[0, 10],
                         npts=npts)

        CouplingLP01.append( PointFieldCoupling(LP01.Fourier, Scat.Sample.Total))
        CouplingLP11.append( PointFieldCoupling(LP11.Fourier, Scat.Sample.Total))
        bar.next()




fig = plt.figure()
ax0 = fig.add_subplot(121)
ax1 = fig.add_subplot(122)

ax0.plot(DiameterList*1e6, CouplingLP01, 'C0', label=r'LP$_{01}$')
ax0.plot(DiameterList*1e6, CouplingLP11, 'C1', label=r'LP$_{11+}$')
ax0.set_title('Mode coupling vs. Scatterer diameter')
ax0.legend()
ax0.set_xlabel(r'Scatterer diameter [$\mu m$]')
ax0.set_ylabel('Modal Coupling')
ax0.grid()

ax1.plot(DiameterList*1e6, np.array(CouplingLP11)/np.array(CouplingLP01))
ax1.set_ylabel('Couplign ratio')
ax1.grid()

plt.show()










# -
