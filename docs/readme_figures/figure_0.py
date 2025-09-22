import numpy as np
from TypedUnit import ureg
import os

dir_path = os.path.dirname(os.path.realpath(__file__))

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup

source = Gaussian(
    wavelength=np.linspace(400, 1000, 500) * ureg.nanometer,
    polarization=0 * ureg.degree,
    optical_power=1e-3 * ureg.watt,
    NA=0.2 * ureg.AU,
)

scatterer = Sphere(
    diameter=[200] * ureg.nanometer,
    property=[4] * ureg.RIU,
    medium_property=1 * ureg.RIU,
    source=source,
)

experiment = Setup(scatterer=scatterer, source=source)
df = experiment.get("Qsca")
df.plot(
    x="source:wavelength", show=False, save_as=f"{dir_path}/../images/resonances.png"
)
