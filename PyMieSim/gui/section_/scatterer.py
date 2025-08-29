from PyMieSim.gui.section_.base import Section, BaseSubSection


class SphereSection(BaseSubSection):
    """Sphere scatterer subsection."""

    name: str = "Sphere"

    def __init__(self, app, parent_section):
        self.app = app
        self.parent_section = parent_section
        self.name_prefix = "scatterer"

        self.inputs = {
            "diameter": {
                "id": "sphere-diameter",
                "label": f"Diameter [nm]",
                "default": "100:20000:200",
            },
            "property": {
                "id": "sphere-property",
                "label": "Scatterer Property",
                "default": "1.5, 1.6",
            },
            "medium_property": {
                "id": "sphere-medium-property",
                "label": "Medium Property",
                "default": "1.33",
            },
        }


class CylinderSection(BaseSubSection):
    """Cylinder scatterer subsection."""

    name: str = "Cylinder"

    def __init__(self, app, parent_section):
        self.app = app
        self.parent_section = parent_section
        self.name_prefix = "scatterer"

        self.inputs = {
            "diameter": {
                "id": "cylinder-diameter",
                "label": f"Diameter [nm]",
                "default": "100:20000:200",
            },
            "property": {
                "id": "cylinder-property",
                "label": "Scatterer Property",
                "default": "1.5, 1.6",
            },
            "medium_property": {
                "id": "cylinder-medium-property",
                "label": "Medium Property",
                "default": "1.33",
            },
        }


class CoreShellSection(BaseSubSection):
    """Core-shell scatterer subsection."""

    name: str = "Core-Shell"

    def __init__(self, app, parent_section):
        self.app = app
        self.parent_section = parent_section
        self.name_prefix = "scatterer"

        self.inputs = {
            "core_diameter": {
                "id": "coreshell-core-diameter",
                "label": f"Core Diameter [nm]",
                "default": "50:200:200",
            },
            "shell_thickness": {
                "id": "coreshell-shell-thickness",
                "label": f"Shell Thickness [nm]",
                "default": "200",
            },
            "core_property": {
                "id": "coreshell-core-property",
                "label": "Core Property",
                "default": "2.4, 2.5",
            },
            "shell_property": {
                "id": "coreshell-shell-property",
                "label": "Shell Property",
                "default": "1.5, 1.6",
            },
            "medium_property": {
                "id": "coreshell-medium-property",
                "label": "Medium Property",
                "default": "1.33",
            },
        }


class ScattererSection(Section):
    """
    Main scatterer section that switches between sphere and core-shell subsections.
    """

    name = "Scatterer"

    def __init__(self, app):
        """Initialize the ScattererSection with subsections."""
        self.color = "green"

        # Create subsections (pass self as parent_section)
        self.sub_sections = [
            SphereSection(app, self),
            CoreShellSection(app, self),
            CylinderSection(app, self),
        ]

        super().__init__(app)
