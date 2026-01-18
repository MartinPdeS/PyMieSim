from PyMieSim.gui.section_.base import Section, BaseSubSection


class PlaneWaveSection(BaseSubSection):
    """Sphere scatterer subsection."""

    name: str = "PlaneWave"

    def __init__(self, app, parent_section):
        self.app = app
        self.parent_section = parent_section
        self.name_prefix = "source"

        self.inputs = {
            "wavelength": {
                "id": "planewave-wavelength",
                "label": f"Wavelength [nm]",
                "default": "750",
            },
            "amplitude": {
                "id": "planewave-property",
                "label": "amplitude [volt/meter]",
                "default": "10",
            },
            "polarization": {
                "id": "planewave-polarization",
                "label": "Polarization [degrees]",
                "default": "0",
            },
        }


class GaussianBeamSection(BaseSubSection):
    """Sphere scatterer subsection."""

    name: str = "Gaussian"

    def __init__(self, app, parent_section):
        self.app = app
        self.parent_section = parent_section
        self.name_prefix = "source"

        self.inputs = {
            "wavelength": {
                "id": "gaussian-wavelength",
                "label": f"Wavelength [nm]",
                "default": "750",
            },
            "numerical_aperture": {
                "id": "gaussian-numerical_aperture",
                "label": "Numerical Aperture",
                "default": "0.1",
            },
            "optical_power": {
                "id": "gaussian-optical_power",
                "label": "Optical Power [mW]",
                "default": "10",
            },
            "polarization": {
                "id": "gaussian-polarization",
                "label": "Polarization [degrees]",
                "default": "0",
            },
        }


class SourceSection(Section):
    """
    Main scatterer section that switches between sphere and core-shell subsections.
    """

    name = "Source"

    def __init__(self, app):
        """Initialize the SourceSection with subsections."""
        self.color = "blue"

        # Create subsections (pass self as parent_section)
        self.sub_sections = [
            PlaneWaveSection(app, self),
            GaussianBeamSection(app, self),
        ]

        super().__init__(app)
