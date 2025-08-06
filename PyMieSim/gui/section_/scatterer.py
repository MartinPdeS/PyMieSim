import numpy
from PyMieSim.gui.section_.base import Section, length_units
from dash import Input, Output, html, dcc
from PyMieSim.gui.helper import parse_string_to_array_or_float

class Base:
    def update_x_axis_options(self, input_values):
        """Update x-axis options for sphere subsection."""
        self._xaxis_options = []
        self._xaxis_options_length = []

        for key, value in input_values.items():
            if value:  # Only process non-empty values
                try:
                    parsed_value = parse_string_to_array_or_float(value)
                    if isinstance(parsed_value, numpy.ndarray) and parsed_value.size > 1:
                        self._xaxis_options.append(f"{self.name_prefix}:{key}")
                        self._xaxis_options_length.append(len(parsed_value))
                except ValueError:
                    pass  # Ignore invalid inputs

    def create_layout(self):
        """Create the sphere input layout."""
        layout = []
        for input_def in self.inputs.values():
            layout.append(html.Div([
                html.Label(input_def['label'], style={'margin-right': '10px'}),
                dcc.Input(id=input_def['id'], type='text', value=input_def['default'], style={'width': '200px'})
            ], style={'margin-bottom': '10px'}))
        return layout

    def get_callback_inputs(self):
        """Get Input objects for callbacks."""
        return [Input(input_def["id"], "value") for input_def in self.inputs.values()]

    def update_callbacks(self):
        """Set up callbacks specific to sphere subsection."""
        self._setup_data_callbacks()


class SphereSection(Base):
    """Sphere scatterer subsection."""

    def __init__(self, name_prefix, app, parent_section):
        self.name_prefix = name_prefix
        self.app = app
        self.parent_section = parent_section
        self._xaxis_options = []
        self._xaxis_options_length = []

        self.inputs = {
            "diameter": {
                "id": f"{name_prefix}-diameter",
                "label": f"Diameter [{length_units}]",
                "default": "100:20000:200"
            },
            "property": {
                "id": f"{name_prefix}-property",
                "label": "Scatterer Property",
                "default": "1.5, 1.6"
            },
            "medium_property": {
                "id": f"{name_prefix}-medium-property",
                "label": "Medium Property",
                "default": "1.33"
            }
        }

    def _setup_data_callbacks(self):
        """Set up data update callbacks for sphere subsection."""
        sphere_inputs = [Input(f"{self.name_prefix}-dropdown", "value")] + self.get_callback_inputs()

        @self.app.callback(
            Output(f"{self.name_prefix}-sphere-data", "children"),
            sphere_inputs,
            prevent_initial_call=False
        )
        def update_sphere_data(scatterer_type, *input_values):
            """Update data for sphere subsection."""
            if scatterer_type != 'sphere':
                return "Not Active"

            # Create input values dictionary
            input_values_dict = {}
            for i, (key, _) in enumerate(self.inputs.items()):
                if i < len(input_values):
                    input_values_dict[key] = input_values[i]

            print(f"Sphere data update - Inputs: {input_values_dict}")

            # Update x-axis options
            self.update_x_axis_options(input_values_dict)

            # Update parent section data
            self.parent_section.data = {'type': 'sphere', **input_values_dict}
            self.parent_section._xaxis_options = self._xaxis_options
            self.parent_section._xaxis_options_length = self._xaxis_options_length

            return "Sphere Data Updated"

class CylinderSection(Base):
    """Cylinder scatterer subsection."""

    def __init__(self, name_prefix, app, parent_section):
        self.name_prefix = name_prefix
        self.app = app
        self.parent_section = parent_section
        self._xaxis_options = []
        self._xaxis_options_length = []

        self.inputs = {
            "diameter": {
                "id": f"{name_prefix}-diameter",
                "label": f"Diameter [{length_units}]",
                "default": "100:20000:200"
            },
            "property": {
                "id": f"{name_prefix}-property",
                "label": "Scatterer Property",
                "default": "1.5, 1.6"
            },
            "medium_property": {
                "id": f"{name_prefix}-medium-property",
                "label": "Medium Property",
                "default": "1.33"
            }
        }

    def _setup_data_callbacks(self):
        """Set up data update callbacks for cylinder subsection."""
        cylinder_inputs = [Input(f"{self.name_prefix}-dropdown", "value")] + self.get_callback_inputs()

        @self.app.callback(
            Output(f"{self.name_prefix}-cylinder-data", "children"),
            cylinder_inputs,
            prevent_initial_call=False
        )
        def update_cylinder_data(scatterer_type, *input_values):
            """Update data for cylinder subsection."""
            if scatterer_type != 'cylinder':
                return "Not Active"

            # Create input values dictionary
            input_values_dict = {}
            for i, (key, _) in enumerate(self.inputs.items()):
                if i < len(input_values):
                    input_values_dict[key] = input_values[i]

            print(f"Cylinder data update - Inputs: {input_values_dict}")

            # Update x-axis options
            self.update_x_axis_options(input_values_dict)

            # Update parent section data
            self.parent_section.data = {'type': 'sphere', **input_values_dict}
            self.parent_section._xaxis_options = self._xaxis_options
            self.parent_section._xaxis_options_length = self._xaxis_options_length

            return "Sphere Data Updated"


class CoreShellSection(Base):
    """Core-shell scatterer subsection."""

    def __init__(self, name_prefix, app, parent_section):
        self.name_prefix = name_prefix
        self.app = app
        self.parent_section = parent_section
        self._xaxis_options = []
        self._xaxis_options_length = []

        self.inputs = {
            "core_diameter": {
                "id": f"{name_prefix}-core-diameter",
                "label": f"Core Diameter [{length_units}]",
                "default": "50:200:200"
            },
            "shell_thickness": {
                "id": f"{name_prefix}-shell-thickness",
                "label": f"Shell Thickness [{length_units}]",
                "default": "200"
            },
            "core_property": {
                "id": f"{name_prefix}-core-property",
                "label": "Core Property",
                "default": "2.4, 2.5"
            },
            "shell_property": {
                "id": f"{name_prefix}-shell-property",
                "label": "Shell Property",
                "default": "1.5, 1.6"
            },
            "medium_property": {
                "id": f"{name_prefix}-medium-property",
                "label": "Medium Property",
                "default": "1.33"
            }
        }

    def _setup_data_callbacks(self):
        """Set up data update callbacks for core-shell subsection."""
        coreshell_inputs = [Input(f"{self.name_prefix}-dropdown", "value")] + self.get_callback_inputs()

        @self.app.callback(
            Output(f"{self.name_prefix}-coreshell-data", "children"),
            coreshell_inputs,
            prevent_initial_call=False
        )
        def update_coreshell_data(scatterer_type, *input_values):
            """Update data for core-shell subsection."""
            if scatterer_type != 'coreshell':
                return "Not Active"

            # Create input values dictionary
            input_values_dict = {}
            for i, (key, _) in enumerate(self.inputs.items()):
                if i < len(input_values):
                    input_values_dict[key] = input_values[i]

            print(f"Core-shell data update - Inputs: {input_values_dict}")

            # Update x-axis options
            self.update_x_axis_options(input_values_dict)

            # Update parent section data
            self.parent_section.data = {'type': 'coreshell', **input_values_dict}
            self.parent_section._xaxis_options = self._xaxis_options
            self.parent_section._xaxis_options_length = self._xaxis_options_length

            return "Core-shell Data Updated"


class ScattererSection(Section):
    """
    Main scatterer section that switches between sphere and core-shell subsections.
    """

    name = 'scatterer'

    def __init__(self, app):
        """Initialize the ScattererSection with subsections."""
        self.app = app
        self.title = "Scatterer"
        self.color = "green"

        # Create subsections (pass self as parent_section)
        self.sphere_section = SphereSection(self.name, app, self)
        self.coreshell_section = CoreShellSection(self.name, app, self)
        self.cylinder_section = CylinderSection(self.name, app, self)

        # Dropdown configuration
        self.dropdown_options = [
            {'label': 'Sphere', 'value': 'sphere'},
            {'label': 'Core-Shell', 'value': 'coreshell'},
            {'label': 'Cylinder', 'value': 'cylinder'},
        ]
        self.dropdown_id = f'{self.name}-dropdown'
        self.default_value = 'sphere'

        # Current active section
        self.current_section = self.sphere_section
        self.inputs = self.current_section.inputs

        # Initialize x-axis options
        self._xaxis_options = []
        self._xaxis_options_length = []

        # Initialize data
        self._initialize_data()

        # Set up callbacks
        self.update_callbacks()

    def _initialize_data(self):
        """Initialize data with default values."""
        default_inputs = {}
        for key, input_def in self.current_section.inputs.items():
            default_inputs[key] = input_def['default']

        self.data = {'type': self.default_value, **default_inputs}
        self.current_section.update_x_axis_options(default_inputs)
        self._xaxis_options = self.current_section._xaxis_options
        self._xaxis_options_length = self.current_section._xaxis_options_length


    def update_callbacks(self):
        """Set up all callbacks."""
        # Callback to switch between subsections
        @self.app.callback(
            Output(f"{self.name}-input-container", "children"),
            Input(self.dropdown_id, "value"),
            prevent_initial_call=False
        )
        def switch_subsection(scatterer_type):
            """Switch between sphere and core-shell subsections."""
            if scatterer_type == 'coreshell':
                self.current_section = self.coreshell_section
            else:
                self.current_section = self.sphere_section

            self.inputs = self.current_section.inputs
            print(f"Switched to {scatterer_type} subsection")
            return self.current_section.create_layout()

        # Set up callbacks for both subsections
        self.current_section.update_callbacks()
        self.coreshell_section.update_callbacks()

    def create(self):
        """Create the main layout with dropdown and dynamic input container."""
        # Title and dropdown
        title = html.H2(self.title, style={'color': self.color})
        dropdown = html.Div([
            html.Label(f"Select {self.title} Type:"),
            dcc.Dropdown(id=self.dropdown_id, options=self.dropdown_options, value=self.default_value)
        ], style={'margin-bottom': '20px'})

        # Dynamic input container
        input_container = html.Div(
            id=f"{self.name}-input-container",
            children=self.current_section.create_layout()
        )

        # Main container
        container = html.Div(
            [title, dropdown, input_container],
            style={'padding': '10px', 'border': '1px solid black', 'margin': '10px'}
        )

        # Hidden data divs for both subsections
        hidden_sphere = html.Div(id=f"{self.name}-sphere-data", style={"display": "none"})
        hidden_coreshell = html.Div(id=f"{self.name}-coreshell-data", style={"display": "none"})
        hidden_cylinder = html.Div(id=f"{self.name}-cylinder-data", style={"display": "none"})

        return container, hidden_sphere, hidden_coreshell