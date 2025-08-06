from dash import Input

from PyMieSim.gui.section_.base import Section, length_units, power_units, angle_units

class SourceSection(Section):
    """
    A class to manage the Source section of the optical simulation interface.

    This class handles the inputs and callbacks related to the light source, such as its type,
    wavelength, optical power, numerical aperture, and polarization.

    Parameters
    ----------
    app : Dash
        The Dash application instance.
    """

    name = 'source'

    def __init__(self, app):
        """
        Initialize the SourceSection class.

        Parameters
        ----------
        app : Dash
            The Dash application instance.
        """
        self.app = app
        self.title = "Source"
        self.color = "blue"
        self.dropdown_options = [
            {'label': 'Laser', 'value': 'laser'}
        ]
        self.dropdown_id = f'{self.name}-dropdown'
        self.default_value = 'laser'

        self.inputs = {
            "wavelength": {
                "id": f"{self.name}-wavelength",
                "label": f"Wavelength [{length_units}]",
                "default": "650"
            },
            "optical_power": {
                "id": f"{self.name}-optical_power",
                "label": f"Power [{power_units}]",
                "default": "5"
            },
            "NA": {
                "id": f"{self.name}-NA",
                "label": "Numerical aperture [NA]",
                "default": "0.2"
            },
            "polarization": {
                "id": f"{self.name}-polarization",
                "label": f"Polarization [{angle_units}]",
                "default": "0"
            }
        }

        # Initialize the x-axis options attribute
        self._xaxis_options = []
        self._xaxis_options_length = []

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [
            Input(input_data["id"], "value") for input_data in self.inputs.values()
        ]

        self.update_callbacks()
