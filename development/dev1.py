
import numpy as np

from PyMieSim.experiment.scatterer import Sphere
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, watt, AU, RIU

source = Gaussian(
    wavelength=488 * nanometer,
    polarization=0 * degree,
    optical_power=1e-3 * watt,
    NA=0.2 * AU
)

scatterer = Sphere(
    diameter=np.linspace(30, 300, 200) * nanometer,
    medium_property=[1.33, 1.3299] * RIU,
    property=[1.42] * RIU,
    source=source
)

detector = Photodiode(
    NA=0.2 * AU,
    phi_offset=[0.0] * degree,
    gamma_offset=0.0 * degree,
    sampling=600 * AU,
    polarization_filter=None
)


experiment = Setup(scatterer=scatterer, source=source, detector=detector)


dataframe = experiment.get('coupling', scale_unit=True)


dataframe.plot_data(x='scatterer:diameter', log_scale_y=False)
