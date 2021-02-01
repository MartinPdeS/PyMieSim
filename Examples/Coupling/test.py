import numpy as np
import matplotlib.pyplot as plt
from PyMieCoupling.classes.Detector import LPmode, Photodiode
from PyMieCoupling.utils import Source
from PyMieCoupling.classes.Sets import ScattererSet, ExperimentalSet
from PyMieCoupling.classes.Scattering import Scatterer

LightSource = Source(Wavelength = 450e-9,
                     Polarization = 0,
                     Power = 1,
                     Radius = 1)

Detector0 = Photodiode(NA                = 0.1,
                       Sampling          = 801,
                       GammaOffset       = 0,
                       PhiOffset         = 180,
                       CouplingMode      = 'Centered')

Detector1 = LPmode(NA                = 0.1,
                   Sampling          = 401,
                   GammaOffset       = 0,
                   PhiOffset         = 0,
                   Mode              = (1,1),
                   CouplingMode      = 'Centered')



Scat = Scatterer(Index=1.4, Source = LightSource, Diameter=1000e-9)

Para, Perp = Detector1.Footprint(Scatterer = Scat, Num=100)
Perp = Perp*0
n=100
FourierPara = np.fft.ifft2(Para,)
FourierPara = np.fft.fftshift(FourierPara).__abs__()**2

FourierPerp = np.fft.ifft2(Perp,)
FourierPerp = np.fft.fftshift(FourierPerp).__abs__()**2


Fourier = FourierPerp + FourierPara

fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(121)
ax1 = fig.add_subplot(122)
im = ax.pcolormesh(Para.__abs__()**2 + 0*Perp.__abs__()**2)
im1 = ax1.pcolormesh(np.log(FourierPara))
fig.colorbar(im, ax = ax)
fig.colorbar(im1, ax = ax1)
plt.show()
