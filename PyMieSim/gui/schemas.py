"""Declarative field schemas for the experiment dashboard."""

from dataclasses import dataclass
from typing import Any, Optional

from PyMieSim.units import ureg


@dataclass(frozen=True)
class FieldSpec:
    """Describe one configurable experiment field rendered by the dashboard."""

    name: str
    label: str
    kind: str
    default: str
    unit: Optional[Any] = None
    optional: bool = False
    help_text: str = ""
    placeholder: str = ""


SOURCE_FIELDS = {
    "GaussianSet": (
        FieldSpec("wavelength", "Wavelength", "quantity", "600:1000:150", ureg.nanometer, help_text="Use a single value, comma-separated values, or start:end:count."),
        FieldSpec("polarization", "Polarization Angles", "polarization", "0,90", ureg.degree, help_text="Angles in degrees used to build a PolarizationSet."),
        FieldSpec("optical_power", "Optical Power", "quantity", "1e-3", ureg.watt),
        FieldSpec("numerical_aperture", "Numerical Aperture", "numeric", "0.2"),
    ),
    "PlaneWaveSet": (
        FieldSpec("wavelength", "Wavelength", "quantity", "400:1000:350", ureg.nanometer),
        FieldSpec("polarization", "Polarization Angles", "polarization", "0", ureg.degree),
        FieldSpec("amplitude", "Amplitude", "quantity", "1", ureg.volt / ureg.meter),
    ),
}

SCATTERER_FIELDS = {
    "SphereSet": (
        FieldSpec("diameter", "Diameter", "quantity", "100,200", ureg.nanometer),
        FieldSpec("material", "Material", "material", "4", help_text="Use refractive indices, complex values, or names such as fused_silica or silver."),
        FieldSpec("medium", "Medium", "medium", "1.0", help_text="Use refractive indices or named media such as water."),
    ),
    "InfiniteCylinderSet": (
        FieldSpec("diameter", "Diameter", "quantity", "400:1400:3", ureg.nanometer),
        FieldSpec("material", "Material", "material", "1.4"),
        FieldSpec("medium", "Medium", "medium", "1.0"),
    ),
    "CoreShellSet": (
        FieldSpec("core_diameter", "Core Diameter", "quantity", "800:1000:3", ureg.nanometer),
        FieldSpec("shell_thickness", "Shell Thickness", "quantity", "300", ureg.nanometer),
        FieldSpec("core_material", "Core Material", "material", "1.4"),
        FieldSpec("shell_material", "Shell Material", "material", "1.5"),
        FieldSpec("medium", "Medium", "medium", "1.0"),
    ),
}

DETECTOR_FIELDS = {
    "None": (),
    "PhotodiodeSet": (
        FieldSpec("numerical_aperture", "Numerical Aperture", "numeric", "0.2"),
        FieldSpec("gamma_offset", "Gamma Offset", "quantity", "0", ureg.degree),
        FieldSpec("phi_offset", "Phi Offset", "quantity", "0", ureg.degree),
        FieldSpec("sampling", "Sampling", "integer", "200"),
        FieldSpec("polarization_filter", "Polarization Filter", "quantity", "", ureg.degree, optional=True),
        FieldSpec("medium", "Detector Medium", "medium", "", optional=True),
    ),
    "CoherentModeSet": (
        FieldSpec("mode_number", "Mode Number", "mode", "LP01", help_text="Use a single mode label or comma-separated labels such as LP01,HG11."),
        FieldSpec("numerical_aperture", "Numerical Aperture", "numeric", "0.1"),
        FieldSpec("gamma_offset", "Gamma Offset", "quantity", "0", ureg.degree),
        FieldSpec("phi_offset", "Phi Offset", "quantity", "0", ureg.degree),
        FieldSpec("rotation", "Rotation", "quantity", "0", ureg.degree),
        FieldSpec("sampling", "Sampling", "integer", "200"),
        FieldSpec("polarization_filter", "Polarization Filter", "quantity", "", ureg.degree, optional=True),
        FieldSpec("medium", "Detector Medium", "medium", "", optional=True),
    ),
}

SECTION_FIELDS = {
    "source": SOURCE_FIELDS,
    "scatterer": SCATTERER_FIELDS,
    "detector": DETECTOR_FIELDS,
}


SINGLE_SOURCE_FIELDS = {
    "Gaussian": (
        FieldSpec("wavelength", "Wavelength", "quantity", "632.8", ureg.nanometer),
        FieldSpec("polarization", "Polarization", "quantity", "0", ureg.degree),
        FieldSpec("optical_power", "Optical Power", "quantity", "1e-3", ureg.watt),
        FieldSpec("numerical_aperture", "Numerical Aperture", "numeric", "0.2"),
    ),
    "PlaneWave": (
        FieldSpec("wavelength", "Wavelength", "quantity", "632.8", ureg.nanometer),
        FieldSpec("polarization", "Polarization", "quantity", "0", ureg.degree),
        FieldSpec("amplitude", "Amplitude", "quantity", "1", ureg.volt / ureg.meter),
    ),
}

SINGLE_SCATTERER_FIELDS = {
    "Sphere": (
        FieldSpec("diameter", "Diameter", "quantity", "200", ureg.nanometer),
        FieldSpec("material", "Material", "material", "1.4"),
        FieldSpec("medium", "Medium", "medium", "1.0"),
    ),
    "InfiniteCylinder": (
        FieldSpec("diameter", "Diameter", "quantity", "200", ureg.nanometer),
        FieldSpec("material", "Material", "material", "1.4"),
        FieldSpec("medium", "Medium", "medium", "1.0"),
    ),
    "CoreShell": (
        FieldSpec("core_diameter", "Core Diameter", "quantity", "120", ureg.nanometer),
        FieldSpec("shell_thickness", "Shell Thickness", "quantity", "40", ureg.nanometer),
        FieldSpec("core_material", "Core Material", "material", "1.4"),
        FieldSpec("shell_material", "Shell Material", "material", "1.5"),
        FieldSpec("medium", "Medium", "medium", "1.0"),
    ),
}
