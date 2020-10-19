import matplotlib.pyplot as plt
import numpy as np

from miecoupling.src.classes.Fiber import fiber
from miecoupling.src.classes.Modes import mode
from miecoupling.src.classes.Scattering import Scatterer

from miecoupling.src.functions.couplings import PointFieldCoupling, MeanFieldCoupling



npts = 101

Fiber = fiber(core_radius = 4.2e-6,
              core_index  = 1.4456,
              clad_radius = 20.5e-6,
              clad_index  = 1.4444)

Mode = mode(fiber      = Fiber,
            LPmode     = (1, 1),
            wavelength = 400e-9,
            npts       = npts,
            )

Mode.magnificate(magnification=2.)

Mode.PlotFields()

Scat = Scatterer(diameter      = 500e-9,
                 wavelength    = 400e-9,
                 index         = 1.4,
                 npts          = npts,
                 ThetaBound    = [-180,180],
                 ThetaOffset   = 0,
                 PhiBound      = [-180,180],
                 PhiOffset     = 0)

Scat.GenField(PolarizationAngle=0)

Scat.PlotS1S2()

Scat.Field.PlotStokes(RectangleTheta=[-20,20], RectanglePhi=[-20,20])


#Scat.Sample.PlotStokes()















# -
