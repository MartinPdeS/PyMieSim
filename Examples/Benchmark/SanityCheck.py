import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from PyMieCoupling.classes.Detector import Photodiode
from PyMieCoupling.classes.Scattering import Scatterer, ScattererSet
from PyMieCoupling.classes.Fields import Source
from scipy import interpolate
from scipy.integrate import dblquad



LightSource = Source(Wavelength   = 450e-9,
                     Polarization = 0)


Photodiode0 = Photodiode(NA                = 0.1,
                         Source            = LightSource,
                         Npts              = 10,
                         ThetaOffset       = 0,
                         PhiOffset         = 0)

DiameterList  = np.linspace(10e-9, 5000e-9, 3)


Set1 = ScattererSet(DiameterList  = DiameterList,
                    RIList        = [1.4],
                    Detector      = Photodiode0,
                    Source        = LightSource,
                    Mode          = 'Centered'
                    )


DF1 = Set1.GetFrame(Filter=0)
DF1.Plot('Coupling')


Theoretical_coupling = []
for diameter in DiameterList:

    Scat = Scatterer(Diameter      = diameter,
                     Source        = LightSource,
                     Index         = 1.4,
                     Meshes        = Photodiode0.Meshes)

    S1 = interpolate.interp1d(Scat.Meshes.Phi.Vector.Radian, Scat.S1S2.S1S2[0])
    S2 = interpolate.interp1d(Scat.Meshes.Phi.Vector.Radian, Scat.S1S2.S1S2[1])



    def integrand(theta, phi):
        temp0 = np.abs(S2(phi))**2 * np.cos(theta)**2

        temp1 = np.abs(S1(phi))**2 * np.sin(theta)**2

        return (temp0 + temp1)/Scat.Source.k**2 * np.abs(np.sin(phi))


    ans, err = dblquad(integrand,
                       *Scat.Meshes.Phi.Boundary.Radian,                  #phi
                       *Scat.Meshes.Theta.Boundary.Radian)                #theta

    Theoretical_coupling.append(ans)










data0 = DF1.Coupling.to_numpy()
data0 /= np.max(data0)


data1 = Theoretical_coupling
data1 /= np.max(data1)
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
