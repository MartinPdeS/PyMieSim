from dash import Input
from PyMieSim.gui.section_.base import Section, angle_units


class DetectorSection(Section):
    """
    A class to manage the Detector section of the optical simulation interface.

    This class handles the inputs and callbacks related to the detector, such as numerical aperture,
    offsets, sampling, and polarization filter.

    Parameters
    ----------
    app : Dash
        The Dash application instance.
    """

    name = 'detector'

    def __init__(self, app):
        """
        Initialize the DetectorSection class.

        Parameters
        ----------
        app : Dash
            The Dash application instance.
        """
        self.app = app
        self.title = "Detector"
        self.color = "red"
        self.dropdown_options = [
            {'label': 'Photodiode', 'value': 'photodiode'},
        ]
        self.dropdown_id = 'detector-dropdown'
        self.default_value = 'photodiode'

        self.inputs = {
            "NA": {
                "id": f"{self.name}-na",
                "label": "Numerical aperture",
                "default": "0.9"
            },
            "phi_offset": {
                "id": f"{self.name}-phi-offset",
                "label": f"Phi offset [{angle_units}]",
                "default": "0"
            },
            "gamma_offset": {
                "id": f"{self.name}-gamma-offset",
                "label": f"Gamma offset [{angle_units}]",
                "default": "0"
            },
            "sampling": {
                "id": f"{self.name}-sampling",
                "label": "Sampling",
                "default": "1000"
            },
            "polarization_filter": {
                "id": f"{self.name}-polarization_filter",
                "label": f"Polarization Filter [{angle_units}]",
                "default": "None"
            }
        }

        # Initialize x-axis options
        self._xaxis_options = []
        self._xaxis_options_length = []

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [
            Input(input_data["id"], "value") for input_data in self.inputs.values()
        ]

        self.update_callbacks()
