try:
    from ._version import version as __version__

except ImportError:
    __version__ = "0.0.0"

import PyMieSim.units as _
import PyMieSim.material as _
import PyMieSim.polarization as _

from .api import Simulation
from .measures import (
    Cabs, Cback, Cext, Cforward, Cpr, Cratio, Csca,
    Measure, MeasureLike, MeasureName, Qabs, Qback, Qext, Qforward, Qpr, Qratio, Qsca,
    coupling, cross_section, g, g_with_farfields, size_parameter,
)
from .results import ExperimentResult, SimulationResult, SimulationResults
from .coordinates import Cartesian, Spherical, VectorField
from .experiment.detector_set import CoherentModeSet, PhotodiodeSet
from .experiment.material_set import MaterialSet, MediumSet
from .experiment.polarization_set import PolarizationSet
from .experiment.scatterer_set import CoreShellSet, InfiniteCylinderSet, SphereSet
from .experiment.source_set import GaussianSet, PlaneWaveSet
from .experiment.setup import Setup as Experiment
from .material import (
    ConstantMaterial,
    ConstantMedium,
    print_available,
    SellmeierMaterial,
    SellmeierMedium,
    TabulatedMaterial,
    TabulatedMedium,
)
from .materials import (
    MaterialInfo,
    available_materials,
    load_material,
    load_tabulated,
    material_info,
    validate_refractive_indices,
    validate_material,
    validate_tabulated_data,
    validate_wavelength,
)
from .mesh import FibonacciMesh, FullMesh
from .labeled_array import LabeledArray
from .polarization import LeftCircular, PolarizationState, RightCircular
from .single import Setup
from .single.detector import CoherentMode, IntegratingSphere, Photodiode
from .single.scatterer import CoreShell, InfiniteCylinder, Sphere
from .single.source import Gaussian, PlaneWave
from .units import ureg

__all__ = [
    "Cartesian",
    "CoherentMode",
    "ConstantMaterial",
    "ConstantMedium",
    "Cabs",
    "Cback",
    "Cext",
    "Cforward",
    "Cpr",
    "Cratio",
    "Csca",
    "CoreShell",
    "CoreShellSet",
    "CoherentModeSet",
    "Experiment",
    "FibonacciMesh",
    "FullMesh",
    "Gaussian",
    "InfiniteCylinder",
    "IntegratingSphere",
    "LeftCircular",
    "LabeledArray",
    "MaterialSet",
    "MaterialInfo",
    "Measure",
    "MeasureLike",
    "MeasureName",
    "MediumSet",
    "Photodiode",
    "PhotodiodeSet",
    "PlaneWave",
    "PolarizationState",
    "PolarizationSet",
    "Qabs",
    "Qback",
    "Qext",
    "Qforward",
    "Qpr",
    "Qratio",
    "Qsca",
    "RightCircular",
    "SellmeierMaterial",
    "SellmeierMedium",
    "Setup",
    "Simulation",
    "Sphere",
    "SphereSet",
    "SimulationResult",
    "SimulationResults",
    "ExperimentResult",
    "Spherical",
    "TabulatedMaterial",
    "TabulatedMedium",
    "VectorField",
    "coupling",
    "cross_section",
    "g",
    "g_with_farfields",
    "size_parameter",
    "GaussianSet",
    "InfiniteCylinderSet",
    "PlaneWaveSet",
    "print_available",
    "available_materials",
    "load_material",
    "load_tabulated",
    "material_info",
    "validate_refractive_indices",
    "validate_material",
    "validate_tabulated_data",
    "validate_wavelength",
    "__version__",
    "ureg",
]
