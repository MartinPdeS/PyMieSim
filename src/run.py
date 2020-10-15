import matplotlib.pyplot as plt
import numpy as np

from classes.Fiber import fiber
from classes.Modes import mode
from classes.Scattering import Scatterer

from functions.couplings import PointFieldCoupling, MeanFieldCoupling

npts=101

Fiber = fiber(4.2e-6,
              1.4456,
              20.5e-6,
              1.4444)

Mode = mode(fiber=Fiber,
            LPmode=(2, 1),
            wavelength=400e-9,
            npts=npts,
            )

Mode.magnificate(magnification=2.)

# Mode.PlotFields()

Scat = Scatterer(diameter=500e-10,
                 wavelength=400e-9,
                 index=1.4,
                 npts=npts)

Scat.GenField(PolarizationAngle=0)

Scat.PlotFields()

Scat.SampleField(ThetaBound=[0,360],
                 PhiBound=[0, 180],
                 npts=npts)


Scat.Sample.PlotStokes()

res = MeanFieldCoupling(Mode.Fourier, Scat.Total)
