import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from PyMieSim.Detector import Photodiode
from PyMieSim.Scatterer import Sphere
from PyMieSim.Sets import ScattererSet, ExperimentalSet
from PyMieSim.Source import PlaneWave
from scipy import interpolate
from scipy.integrate import dblquad
from PyMieSim.LMT.Sphere import S1S2 as GetS1S2


LightSource = PlaneWave(Wavelength = 450e-9, Polarization = 0)


Photodiode0 = Photodiode(NA                = 1.,
                         Sampling          = 101,
                         GammaOffset       = 0,
                         PhiOffset         = 0,
                         CouplingMode      = 'Centered')

DiameterList  = np.linspace(100e-9, 2000e-9, 100)

Set1 = ScattererSet(DiameterList  = DiameterList,
                    RIList        = [1.4],
                    Source        = LightSource)


Experience = ExperimentalSet(ScattererSet = Set1, Detectors=Photodiode0)

DF1 = Experience.DataFrame

Theoretical_coupling = []


def integrand(theta, phi):

    temp0 = np.abs(S2(phi))**2 * np.cos(theta)**2

    temp1 = np.abs(S2(phi))**2 * np.sin(theta)**2

    return (temp0 + temp1) * np.sin(phi) / LightSource.k**2


for nD, diameter in enumerate(DiameterList):
    print(f"PROGRESSION: {nD}/{len(DiameterList)}")

    s1, s2 = GetS1S2(1.4,
                     diameter,
                     LightSource.Wavelength,
                     1.0,
                     Photodiode0.Mesh.Phi.Radian-np.pi/2)

    S1 = interpolate.interp1d(Photodiode0.Mesh.Phi.Radian, s1)
    S2 = interpolate.interp1d(Photodiode0.Mesh.Phi.Radian, s2)

    ans, err = dblquad(integrand,
                       Photodiode0.Mesh.Phi.Radian.max(),
                       Photodiode0.Mesh.Phi.Radian.min(),
                       0,
                       2*np.pi,
                       epsabs=1e-8)

    Theoretical_coupling.append(ans)



data0 = DF1.Coupling.to_numpy()
data0 /= np.mean(data0)


data1 = Theoretical_coupling
data1 /= np.mean(data1)
text1 = r"$\iint_\Omega \frac{ |S_2|^2  \cos{\theta}^2 + |S1|^2  \sin{\theta}^2 }{k^2} \sin{\phi} d\phi d\theta$"

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
