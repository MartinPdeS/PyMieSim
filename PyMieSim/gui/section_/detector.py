from PyMieSim.gui.section_.base import Section, BaseSubSection


class PhotodiodeSection(BaseSubSection):
    """Photodiode scatterer subsection."""

    name: str = "photodiode"

    def __init__(self, app, parent_section):
        self.app = app
        self.parent_section = parent_section
        self.name_prefix = "detector"

        self.inputs = {
            "numerical_aperture": {
                "id": f"photodiode-na",
                "label": "Numerical aperture",
                "default": "0.2",
            },
            "phi_offset": {
                "id": f"photodiode-phi-offset",
                "label": f"Phi offset [degrees]",
                "default": "0",
            },
            "gamma_offset": {
                "id": f"photodiode-gamma-offset",
                "label": f"Gamma offset [degrees]",
                "default": "0",
            },
            "sampling": {
                "id": f"photodiode-sampling",
                "label": "Sampling",
                "default": "1000",
            },
            "polarization_filter": {
                "id": f"photodiode-polarization_filter",
                "label": f"Polarization Filter [degrees]",
                "default": "None",
            },
        }


class CoherentModeSection(BaseSubSection):
    """Coherent mode detector subsection."""

    name: str = "coherentmode"

    def __init__(self, app, parent_section):
        self.app = app
        self.parent_section = parent_section
        self.name_prefix = "detector"

        self.inputs = {
            "mode_number": {
                "id": f"coherentmode-field",
                "label": "Mode field",
                "default": "LP01, LP02",
            },
            "numerical_aperture": {
                "id": f"coherentmode-numerical-aperture",
                "label": "Numerical aperture",
                "default": "0.2",
            },
            "phi_offset": {
                "id": f"coherentmode-phi-offset",
                "label": f"Phi offset [degrees]",
                "default": "0",
            },
            "gamma_offset": {
                "id": f"coherentmode-gamma-offset",
                "label": f"Gamma offset [degrees]",
                "default": "0",
            },
            "rotation": {
                "id": f"coherentmode-rotation",
                "label": f"Rotation [degrees]",
                "default": "0",
            },
            "sampling": {
                "id": f"coherentmode-sampling",
                "label": "Sampling",
                "default": "1000",
            },
            "polarization_filter": {
                "id": f"coherentmode-polarization_filter",
                "label": f"Polarization Filter [degrees]",
                "default": "None",
            },
        }


class DetectorSection(Section):
    """
    Main scatterer section that switches between sphere and core-shell subsections.
    """

    name = "Detector"

    def __init__(self, app):
        """Initialize the DetectorSection with subsections."""
        self.color = "blue"

        # Create subsections (pass self as parent_section)
        self.sub_sections = [
            PhotodiodeSection(app, self),
            CoherentModeSection(app, self),
        ]

        super().__init__(app)
