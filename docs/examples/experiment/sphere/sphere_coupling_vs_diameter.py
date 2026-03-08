"""
Sphere: Coupling vs diameter
============================

"""

# %%
# Importing the package dependencies: numpy, PyMieSim
import numpy
from PyMieSim.units import ureg

from PyMieSim.single.material import Material, Medium
from PyMieSim.experiment.detector_set import CoherentModeSet
from PyMieSim.experiment.scatterer_set import SphereSet
from PyMieSim.experiment.source_set import GaussianSet
from PyMieSim.experiment.polarization_set import PolarizationSet
from PyMieSim.experiment.material_set import MaterialSet, MediumSet

from PyMieSim.experiment.setup import Setup
import PyMieSim
PyMieSim.debug_mode = True
# from PyOptik import Material

polarization_set = PolarizationSet(
    angles=[90.0] * ureg.degree,
)

source = GaussianSet(
    wavelength=[700] * ureg.nanometer,
    polarization=polarization_set,
    optical_power=[1e-3] * ureg.watt,
    numerical_aperture=[0.1] * ureg.AU,
)


bk7 = Material(
    name="BK7",
    refractive_indices=[1.3, 1.4, 1.5, 1.6, 1.7] * ureg.RIU,
    wavelengths=[400, 500, 600, 700, 800] * ureg.nanometer,
)


bk2 = Material(
    name="BK2",
    refractive_indices=[1.2, 1.4, 1.5, 1.6, 1.7] * ureg.RIU,
    wavelengths=[400, 500, 600, 700, 800] * ureg.nanometer,
)


water = Medium(
    name="Water",
    refractive_indices=[1.1, 1.2] * ureg.RIU,
    wavelengths=[400, 800] * ureg.nanometer,
)


materials = MaterialSet(
    # [1.4, 1.5]
    [bk7, bk2]
)


mediums = MediumSet(
    # [1.4, 1.5]
    [water]
)


scatterer = SphereSet(
    diameter=numpy.linspace(100, 10000, 600) * ureg.nanometer,
    # material=[Material.BK7],
    material=materials,
    # material=materials,
    medium=mediums,
    source=source,
)

detector = CoherentModeSet(
    mode_number=["LP01", "LP11", "LP02"],
    numerical_aperture=[0.2] * ureg.AU,
    rotation=[0] * ureg.degree,
    phi_offset=[0.0] * ureg.degree,
    gamma_offset=[0.0] * ureg.degree,
    sampling=[600],
    mean_coupling=True,
)

experiment = Setup(scatterer=scatterer, source=source, detector=detector)

dataframe = experiment.get("coupling", drop_unique_level=True, scale_unit=True)

dataframe.plot(x="scatterer:diameter")
