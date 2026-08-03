"""
Sphere: Coupling vs diameter
============================

"""
import numpy
from PyMieSim.units import ureg

from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment import Setup


polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[700] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.1],
)

scatterer = SphereSet(
    diameter=numpy.linspace(100, 10000, 600) * ureg.nanometer,
    material=[1.5],
    medium=[1.0],
)

detector = CoherentModeSet(
    mode_number=["LP01", "LP11", "LP02"],
    numerical_aperture=[0.2],
    rotation=[0] * ureg.degree,
    phi_offset=[0.0] * ureg.degree,
    gamma_offset=[0.0] * ureg.degree,
    sampling=[600],
    mean_coupling=True,
)

experiment = Setup(scatterer_set=scatterer, source_set=source, detector_set=detector)

dataframe = experiment.get("coupling", drop_unique_level=True, scale_unit=True)

dataframe.plot(x="scatterer:diameter")
