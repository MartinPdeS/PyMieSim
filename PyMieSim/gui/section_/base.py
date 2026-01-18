from dash import html, dcc, Input, Output
from TypedUnit import ureg
from PyMieSim.gui.helper import parse_string_to_array_or_float
import numpy

length_units = ureg.nanometer
power_units = ureg.milliwatt
angle_units = ureg.degree


class Section:
    """
    A base class to define shared functionality for different sections in the optical simulation interface.

    This class provides methods for creating dropdown menus, input fields, and setting up callbacks for
    interaction with the Dash application.
    """

    @property
    def label(self):
        return self.name.replace("_", "").replace("-", "").lower()

    def __init__(self, app):
        self.app = app

        # Dropdown configuration
        self.dropdown_id = f"{self.label}-dropdown"
        self.dropdown_options = [
            {"label": f"{sub_section.name}", "value": f"{sub_section.label}"}
            for sub_section in self.sub_sections
        ]

        # Current active section
        self.current_section = self.sub_sections[0]
        self.default_value = self.current_section.label

        self.inputs = self.current_section.inputs

        # Initialize data
        self._initialize_data()

        # Set up callbacks
        self.update_callbacks()

    def _initialize_data(self) -> None:
        """Initialize data with default values."""
        default_inputs = {}
        for key, input_def in self.current_section.inputs.items():
            default_inputs[key] = input_def["default"]

        self.data = {"type": self.default_value, **default_inputs}
        self.current_section.update_x_axis_options(default_inputs)

    def create(self):
        """Create the main layout with dropdown and dynamic input container."""
        # Title and dropdown
        title = html.H2(self.name, style={"color": self.color})
        dropdown = html.Div(
            [
                html.Label(f"Select {self.name} Type:"),
                dcc.Dropdown(
                    id=self.dropdown_id,
                    options=self.dropdown_options,
                    value=self.default_value,
                ),
            ],
            style={"margin-bottom": "20px"},
        )

        # Dynamic input container
        input_container = html.Div(
            id=f"{self.label}-input-container",
            children=self.current_section.create_layout(),
        )

        # Main container
        container = html.Div(
            [title, dropdown, input_container],
            style={"padding": "10px", "border": "1px solid black", "margin": "10px"},
        )

        return container, *self._get_hidden_data_divs()

    def update_callbacks(self) -> None:
        """Set up all callbacks."""
        # Callback to switch between subsections
        api_outputs = Output(f"{self.label}-input-container", "children")
        api_inputs = [Input(self.dropdown_id, "value")]

        @self.app.callback(api_outputs, api_inputs, prevent_initial_call=False)
        def switch_subsection(source_type):
            """Switch between plane wave subsections."""
            return self._switch_subsection(source_type)

        # Set up callbacks for both subsections
        [sub_section._setup_data_callbacks() for sub_section in self.sub_sections]

    def _get_hidden_data_divs(self):
        """Get hidden data divs for all subsections."""
        hidden_subsections = [
            html.Div(
                id=f"{self.label}-{sub_section.label}-data", style={"display": "none"}
            )
            for sub_section in self.sub_sections
        ]
        return hidden_subsections

    def _switch_subsection(self, object_type):
        """Switch between plane wave subsections."""
        for sub_section in self.sub_sections:
            if sub_section.label == object_type:
                self.current_section = sub_section
                break

        self.inputs = self.current_section.inputs
        print(f"Switched to {object_type} subsection")
        return self.current_section.create_layout()


class BaseSubSection:

    @property
    def label(self):
        return self.name.replace("_", "").replace("-", "").lower()

    def update_x_axis_options(self, input_values) -> list:
        """Update x-axis options for sphere subsection."""
        self.parent_section._xaxis_options = []

        for key, value in input_values.items():
            if value:  # Only process non-empty values
                try:
                    parsed_value = parse_string_to_array_or_float(value)
                    if (
                        isinstance(parsed_value, numpy.ndarray)
                        and parsed_value.size > 1
                    ):
                        self.parent_section._xaxis_options.append(f"{self.name}:{key}")
                except ValueError:
                    pass  # Ignore invalid inputs

        self.parent_section._xaxis_options_length = [
            len(x) for x in self.parent_section._xaxis_options
        ]
        self.parent_section.data = {"type": self.label, **input_values}

    def create_layout(self):
        """Create the sphere input layout."""
        layout = []
        for input_def in self.inputs.values():
            layout.append(
                html.Div(
                    [
                        html.Label(input_def["label"], style={"margin-right": "10px"}),
                        dcc.Input(
                            id=input_def["id"],
                            type="text",
                            value=input_def["default"],
                            style={"width": "200px"},
                        ),
                    ],
                    style={"margin-bottom": "10px"},
                )
            )
        return layout

    def _setup_data_callbacks(self) -> None:
        """Set up data update callbacks for cylinder subsection."""
        api_inputs = [Input(f"{self.parent_section.label}-dropdown", "value")]
        api_inputs.extend(
            [Input(input_def["id"], "value") for input_def in self.inputs.values()]
        )

        api_outputs = Output(
            f"{self.parent_section.label}-{self.label}-data", "children"
        )

        @self.app.callback(api_outputs, api_inputs, prevent_initial_call=False)
        def update_scatterer_data(scatterer_type, *input_values):
            """Update data for scatterer subsection."""
            if scatterer_type != self.label:
                return "Not Active"

            # Create input values dictionary
            input_values_dict = {}
            for i, (key, _) in enumerate(self.inputs.items()):
                if i < len(input_values):
                    input_values_dict[key] = input_values[i]

            # Update x-axis options
            self.update_x_axis_options(input_values_dict)
