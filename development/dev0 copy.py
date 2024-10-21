import numpy as np
import matplotlib.pyplot as plt
import PyMieScatt

from PyMieSim.experiment.source import Gaussian, PlaneWave
from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment import Setup
from PyMieSim.polarization import Linear
from PyMieSim.units import nanometer, degree, watt, AU, RIU
from PyOptik import Material

def get_m(material, wl_range):
    n_list = list(material.n_values)
    k_list = list(material.k_values)
    wl_mat_um = list(material.wavelength)
    wl_mat = [wl_mat_um[i]*1e3 for i in range(len(wl_mat_um))]

    # Linear interpolation
    n_interp = np.interp(wl_range, wl_mat, n_list)
    k_interp = np.interp(wl_range, wl_mat, k_list)

    # Create complex index array from interpolated values
    m = np.array([complex(n_interp[i], k_interp[i]) for i in range(len(n_interp))])

    return m

wl_range = np.linspace(200, 1200, 100)  # Wavelength range in nanometers

material = Material.gold
r = 25   # particle radius in nanometers
d = 2*r  # diameter
n_medium = 1.5  # medium index

#%% PyMieSim
source = PlaneWave(
    wavelength=wl_range * nanometer,
    polarization=0 * degree,
    amplitude=0 * watt,
)

scatterer = Sphere(
    diameter=d * nanometer,
    property=material,
    medium_property=n_medium * RIU,
    source=source
)
experiment = Setup(
    scatterer=scatterer,
    source=source,
)

qsca_pymiesim = experiment.get('Qsca', add_units=False).squeeze().tolist()
qabs_pymiesim = experiment.get('Qabs', add_units=False).squeeze().tolist()
qext_pymiesim = experiment.get('Qext', add_units=False).squeeze().tolist()

#%% PyMieScatt
m = get_m(material, wl_range)  # Interpolated complex indices array of the particle material

_, qext_pymiescatt, qsca_pymiescatt, qabs_pymiescatt, *_  = PyMieScatt.MieQ_withWavelengthRange(m, d, n_medium, wl_range)

#%% Plotting
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(1, 1, 1)
ax.set_title("Scattering efficiencies (miepython vs PyMieScatt vs PyMieSim)")

ax.plot(wl_range, qext_pymiesim,   ls='-', lw=0.5, marker='o', markersize=4, c='tab:green', label=r"$Q_{ext}^{PyMieSim}$")
ax.plot(wl_range, qext_pymiescatt, ls='-', lw=0.5, marker='o', markersize=4, c='tab:orange', label=r"$Q_{ext}^{PyMieScatt}$")
ax.plot(wl_range, qabs_pymiesim,   ls='-', lw=0.5, marker='^', markersize=4, c='tab:green', label=r"$Q_{abs}^{PyMieSim}$")
ax.plot(wl_range, qabs_pymiescatt, ls='-', lw=0.5, marker='^', markersize=4, c='tab:orange', label=r"$Q_{abs}^{PyMieScatt}$")
ax.plot(wl_range, qsca_pymiesim,   ls='-', lw=0.5, marker='s', markersize=4, c='tab:green', label=r"$Q_{sca}^{PyMieSim}$")
ax.plot(wl_range, qsca_pymiescatt, ls='-', lw=0.5, marker='s', markersize=4, c='tab:orange', label=r"$Q_{sca}^{PyMieScatt}$")

plt.xlabel("Wavelength (nm)")
plt.ylabel("Efficiency")
plt.legend()
plt.tight_layout()
plt.show()