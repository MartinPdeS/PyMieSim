from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup
from PyMieSim.units import ureg
import numpy as np

polarization = PolarizationSet(
    angles=[0] * ureg.degree
)

source = GaussianSet(
    wavelength=np.linspace(400, 2000, 500) * ureg.nanometer,
    polarization=polarization,
    optical_power=1e-3 * ureg.watt,
    numerical_aperture=0.2 * ureg.AU,
)

scatterer = SphereSet(
    diameter=[200, 300] * ureg.nanometer,
    material=[4 + 1j] * ureg.RIU,
    medium=[1] * ureg.RIU,
)

experiment = Setup(
    scatterer_set=scatterer,
    source_set=source
)
df = experiment.get("Qsca", "Qext")
df.plot(x="source:wavelength")