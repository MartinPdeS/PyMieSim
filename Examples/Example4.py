
"""
_________________________________________________________
Optimisation code for RI dependence.
_________________________________________________________
"""

import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from PyMieCoupling.classes.Detector import Detector
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.functions.couplings import PointFieldCoupling

npts=101



Detector = Detector(size       = 50e-6,
                    wavelength = 400e-9,
                    npts       = npts)



def IterDiameter(DiameterList, index):

    Coupling = []

    for Diameter in tqdm(DiameterList, total = len(DiameterList), desc ="Progress:"):

        Scat = Scatterer(diameter    = Diameter,
                         wavelength  = 400e-9,
                         index       = index,
                         npts        = 101,
                         ThetaBound  = Detector.ThetaBound,
                         ThetaOffset = 0,
                         PhiBound    = Detector.PhiBound,
                         PhiOffset   = 0)

        Coupling.append( PointFieldCoupling(Detector = Detector,
                                            Source   = Scat,
                                            Field = 'Parallel') )


    return Coupling




NPTS_I, NPTS_D = 10, 10

IndexList = np.linspace(1.3,2,NPTS_I)

Res = np.empty([NPTS_D, NPTS_I])

for ni, index in enumerate(IndexList):
    Res[:, ni] = IterDiameter(np.linspace(100,3000,NPTS_D) * 1e-9, index=index)




std = np.std(Res, axis=1)


fig, (ax0, ax1) = plt.subplots(1,2, figsize=(10,5))

ax0.set_title('Detec. response vs. scatterer size')

ax1.set_title('Standard deviation Detec. response of RI vs. Scatterer size')

ax0.set_xlabel('Scatterer size')

ax1.set_xlabel('Scatterer size')

ax0.set_ylabel('Detector response [U.A]')

[ ax0.plot( Res[:, ni], label='index: {0:.2f}'. format(index)) for ni, index in enumerate(IndexList)  ]

ax0.grid()

ax0.legend()

ax1.plot(std)

ax1.grid()

plt.show()












# -
