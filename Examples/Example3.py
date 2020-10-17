import matplotlib.pyplot as plt
import numpy as np

from miecoupling.src.classes.Fiber import fiber
from miecoupling.src.classes.Modes import mode
from miecoupling.src.classes.Scattering import Scatterer

from miecoupling.src.functions.couplings import PointFieldCoupling, MeanFieldCoupling



npts=101

Fiber = fiber(4.2e-6,
              1.4456,
              20.5e-6,
              1.4444)

Mode = mode(fiber=Fiber,
            LPmode=(1, 1),
            wavelength=400e-9,
            npts=npts,
            )

Mode.magnificate(magnification=2.)

Mode.PlotFields()

Scat = Scatterer(diameter=1500e-9,
                 wavelength=400e-9,
                 index=1.4,

                 npts=200)



Scat.GenField(PolarizationAngle=0)



Scat.PlotFields()


Scat.SampleField(ThetaBound=[-180,180],
                 PhiBound=[-180, 180],
                 npts=npts)

Scat.Field.PlotStokes(RectangleTheta=[-20,20], RectanglePhi=[-20,20])


#Scat.Sample.PlotStokes()

res = PointFieldCoupling(Mode.Fourier, Scat.Sample.Total)














# -
