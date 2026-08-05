try:
    from ._version import version as __version__

except ImportError:
    __version__ = "0.0.0"

import PyMieSim.units as _
import PyMieSim.material as _
import PyMieSim.polarization as _

from .api import Simulation
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
from .mesh import FibonacciMesh, FullMesh
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
    "MaterialSet",
    "MediumSet",
    "Photodiode",
    "PhotodiodeSet",
    "PlaneWave",
    "PolarizationState",
    "PolarizationSet",
    "RightCircular",
    "SellmeierMaterial",
    "SellmeierMedium",
    "Setup",
    "Simulation",
    "Sphere",
    "SphereSet",
    "Spherical",
    "TabulatedMaterial",
    "TabulatedMedium",
    "VectorField",
    "GaussianSet",
    "InfiniteCylinderSet",
    "PlaneWaveSet",
    "print_available",
    "__version__",
    "ureg",
]
