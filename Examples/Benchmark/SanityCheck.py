import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer
from PyMieCoupling.classes.Sets import ScattererSet, ExperimentalSet
from PyMieCoupling.utils import Source
from scipy import interpolate
from scipy.integrate import dblquad
from PyMieCoupling.cython.S1S2 import GetS1S2


LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)


Photodiode0 = Photodiode(NA                = 0.2,
                         Sampling          = 201,
                         GammaOffset       = 0,
                         PhiOffset         = 0,
                         CouplingMode      = 'Centered')

DiameterList  = np.linspace(100e-9, 5000e-9, 100)

#Photodiode0.Plot()


Set1 = ScattererSet(DiameterList  = DiameterList,
                    RIList        = [1.4],
                    Source        = LightSource)


Experience = ExperimentalSet(ScattererSet = Set1, Detectors=Photodiode0)

DF1 = Experience.DataFrame

Theoretical_coupling = []
for diameter in DiameterList:

    Scat = Scatterer(Diameter      = diameter,
                     Source        = LightSource,
                     Index         = 1.4)

    s1s2 = GetS1S2(Scat.Index, Scat.SizeParam, Photodiode0.Meshes.Phi.Radian)

    S1 = interpolate.interp1d(Scat.S1S2()['Phi'], Scat.S1S2()['S1'])
    S2 = interpolate.interp1d(Scat.S1S2()['Phi'], Scat.S1S2()['S2'])



    def integrand(theta, phi):
        temp0 = np.abs(S2(phi))**2 * np.cos(theta)**2

        temp1 = np.abs(S1(phi))**2 * np.sin(theta)**2

        return (temp0 + temp1) * np.sin(phi) / Scat.Source.k**2


    ans, err = dblquad(integrand,
                       *(0,0.2),                  #phi
                       *(0, np.pi),
                       epsabs=1e-4)                                        #theta

    Theoretical_coupling.append(ans)










data0 = DF1.Coupling.to_numpy()
data0 /= np.mean(data0)


data1 = Theoretical_coupling
data1 /= np.mean(data1)
text1 = r"$\int \int_\Omega \frac{ |S_2|^2  \cos{\theta}^2 + |S1|^2  \sin{\theta}^2 }{k^2} \sin{\phi} d\phi d\theta$"

fig = plt.figure(figsize=(10,5))
ax1 = fig.add_subplot(111)
ax1.plot(DiameterList, data0, 'C0s', label='Simulations')
ax1.plot(DiameterList, data1, "C1", label = text1)
ax1.set_xlabel('Scatterer size [m]')
ax1.set_ylabel('Detector coupling [u.a.]')
ax1.legend()

fig.tight_layout()

ax1.grid()
plt.show()
