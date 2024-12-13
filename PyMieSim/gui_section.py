from dash import html, dcc, State, Input, Output
from PyMieSim.units import nanometer, degree, milliwatt


length_units = nanometer
power_units = milliwatt
angle_units = degree

class Section:
    """
    A base class to define shared functionality for different sections in the optical simulation interface.

    This class provides methods for creating dropdown menus, input fields, and setting up callbacks for
    interaction with the Dash application.
    """

    def create_call_backs(self):
        """
        Create and set up callbacks for the section.

        This method defines the `call_backs` attribute by combining the dropdown and input field callbacks
        and calls the `update_callbacks` method to set up the necessary Dash callbacks.

        Attributes
        ----------
        call_backs : list
            A list of Dash `Input` objects corresponding to the dropdown and input fields of the section.
        """
        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [Input(id, "value") for id in self.input_id]
        self.update_callbacks()

    def create(self):
        """
        Create the layout for the section.

        This method generates the HTML layout for the section, including a title, a dropdown menu to select
        the section type, and input fields for relevant parameters.

        Returns
        -------
        tuple
            A tuple containing:
            - html.Div : The main layout for the section.
            - html.Div : A hidden placeholder for storing section-related data.
        """
        layout = [html.H2(self.title, style={'color': self.color})]

        layout.append(html.Div([
            html.Label(f"Select {self.title} Type:"),
            html.Div([
                dcc.Dropdown(
                    id=self.dropdown_id,
                    options=self.dropdown_options,
                    value=self.default_value
                )
            ])
        ], style={'margin-bottom': '20px'}))

        for input in self.inputs.values():
            layout.append(html.Div([
                html.Label(input['label'], style={'margin-right': '10px'}),
                dcc.Input(
                    id=input['id'],
                    type='text',
                    value=input['default'],
                    style={'width': '200px'}
                )
            ], style={'margin-bottom': '10px'}))

        return (
            html.Div(layout, style={'padding': '10px', 'border': '1px solid black', 'margin': '10px'}),
            html.Div(id=f"{self.name}-dropdown-data", style={"display": "none"})  # Hidden placeholder for data
        )


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

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [
            Input(input_data["id"], "value") for input_data in self.inputs.values()
        ]

        self.update_callbacks()

    def update_callbacks(self):
        """
        Set up Dash callbacks to update the data dictionary based on user input.

        This method listens for changes in the dropdown and input fields and updates
        the `data` attribute with the current values.

        Updates
        -------
        type : str
            The selected source type (e.g., 'Laser').
        wavelength : str
            The wavelength of the light source in nanometers.
        optical_power : str
            The optical power of the light source in milliwatts.
        NA : str
            The numerical aperture of the light source.
        polarization : str
            The polarization of the light source in degrees.
        """
        @self.app.callback(Output(f"{self.dropdown_id}-data", "children"), *self.call_backs)
        def update_data(source_type, wavelength, optical_power, NA, polarization):
            self.data = {
                'type': source_type,
                'wavelength': wavelength,
                'optical_power': optical_power,
                'NA': NA,
                'polarization': polarization
            }
            return "Data Updated"


class ScattererSection(Section):
    """
    A class to manage the Scatterer section of the optical simulation interface.

    This class handles the inputs and callbacks related to the scatterer, such as its type,
    diameter, refractive index property, and medium refractive index.

    Parameters
    ----------
    app : Dash
        The Dash application instance.
    """

    name = 'scatterer'

    def __init__(self, app):
        """
        Initialize the ScattererSection class.

        Parameters
        ----------
        app : Dash
            The Dash application instance.
        """
        self.app = app
        self.title = "Scatterer"
        self.color = "green"
        self.dropdown_options = [
            {'label': 'Sphere', 'value': 'sphere'},
        ]
        self.dropdown_id = 'scatterer-dropdown'
        self.default_value = 'sphere'

        self.inputs = {
            "diameter": {
                "id": f"{self.name}-diameter",
                "label": f"Diameter [{length_units}]",
                "default": "100:20000:200"
            },
            "property": {
                "id": f"{self.name}-property",
                "label": "Scatterer Property",
                "default": "1.5, 1.6"
            },
            "medium_property": {
                "id": f"{self.name}-medium-property",
                "label": "Medium Property",
                "default": "1.33"
            }
        }

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [
            Input(input_data["id"], "value") for input_data in self.inputs.values()
        ]

        self.update_callbacks()

    def update_callbacks(self):
        """
        Set up Dash callbacks to update the data dictionary based on user input.

        This method listens for changes in the dropdown and input fields and updates
        the `data` attribute with the current values.

        Updates
        -------
        type : str
            The selected scatterer type (e.g., 'Sphere').
        diameter : str
            The scatterer diameter value, which can be a range or a single value.
        property : str
            The scatterer refractive index property.
        medium_property : str
            The refractive index of the medium.
        """
        @self.app.callback(Output(f"{self.dropdown_id}-data", "children"), *self.call_backs)
        def update_data(scatterer_type, diameter, property, medium_property):
            self.data = {
                'type': scatterer_type,
                'diameter': diameter,
                'property': property,
                'medium_property': medium_property
            }
            return "Data Updated"


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

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [
            Input(input_data["id"], "value") for input_data in self.inputs.values()
        ]

        self.update_callbacks()

    def update_callbacks(self):
        """
        Set up Dash callbacks to update the data dictionary based on user input.

        This method listens for changes in the dropdown and input fields and updates
        the `data` attribute with the current values.

        Updates
        -------
        type : str
            The selected detector type (e.g., 'Photodiode').
        NA : str
            The numerical aperture value.
        phi_offset : str
            The phi offset value in degrees.
        gamma_offset : str
            The gamma offset value in degrees.
        sampling : str
            The sampling value.
        polarization_filter : str
            The polarization filter value.
        """
        @self.app.callback(Output(f"{self.dropdown_id}-data", "children"), *self.call_backs)
        def update_data(detector_type, NA, phi_offset, gamma_offset, sampling, polarization_filter):
            self.data = {
                'type': detector_type,
                'NA': NA,
                'phi_offset': phi_offset,
                'gamma_offset': gamma_offset,
                'sampling': sampling,
                'polarization_filter': polarization_filter
            }
            return "Data Updated"


class MeasureSection:
    """
    A class to manage the Measure section of the optical simulation interface.

    Parameters
    ----------
    app : Dash
        The Dash application instance.
    scatterer_section : Section
        The Scatterer section instance providing scatterer-related inputs.
    source_section : Section
        The Source section instance providing source-related inputs.
    detector_section : Section
        The Detector section instance providing detector-related inputs.
    """

    def __init__(self, app, scatterer_section, source_section, detector_section):
        self.app = app
        self.scatterer_section = scatterer_section
        self.source_section = source_section
        self.detector_section = detector_section
        self.dropdown_id = "measure-input"
        self.plot_button_id = "generate-plot"
        self.xaxis_input_id = "xaxis-input"
        self.filename_input_id = "filename-input"
        self.save_button_id = "save-data"
        self.plot_ready_store_id = "plot-ready"
        self.data = "Qsca"  # Default measure value

    def get_measure_dropdown(self) -> dcc.Dropdown:
        """
        Creates a dropdown menu for selecting the measure type.

        Returns
        -------
        dcc.Dropdown
            A Dash dropdown component for selecting measures such as Qsca, Qext, etc.
        """
        return dcc.Dropdown(
            id=self.dropdown_id,
            options=[
                {'label': 'Qsca', 'value': 'Qsca'},
                {'label': 'Qext', 'value': 'Qext'},
                {'label': 'Qabs', 'value': 'Qabs'},
                {'label': 'Qpr', 'value': 'Qpr'},
                {'label': 'Csca', 'value': 'Csca'},
                {'label': 'Cext', 'value': 'Cext'},
                {'label': 'Cabs', 'value': 'Cabs'},
                {'label': 'Cpr', 'value': 'Cpr'},
                {'label': 'g', 'value': 'g'},
                {'label': 'Coupling', 'value': 'coupling'}
            ],
            value='Qsca',
            style={'margin-right': '10px', 'margin-bottom': '0px', 'height': '36px', 'width': '200px'}
        )

    def get_xaxis_dropdown(self) -> dcc.Dropdown:
        """
        Creates a dropdown menu for selecting the x-axis parameter.

        Returns
        -------
        dcc.Dropdown
            A Dash dropdown component populated dynamically with x-axis options
            based on the Scatterer, Source, and Detector sections.
        """
        options = []
        options.extend(
            [{"label": "scatterer:" + k, "value": "scatterer:" + k} for k in self.scatterer_section.inputs.keys()]
        )
        options.extend(
            [{"label": "source:" + k, "value": "source:" + k} for k in self.source_section.inputs.keys()]
        )
        options.extend(
            [{"label": "detector:" + k, "value": "detector:" + k} for k in self.detector_section.inputs.keys()]
        )

        default_value = options[0]["value"] if options else None

        return dcc.Dropdown(
            id=self.xaxis_input_id,
            options=options,
            value=default_value,
            style={'margin-right': '10px', 'height': '36px', 'width': '300px'}
        )

    def create(self) -> html.Div:
        """
        Creates the layout for the Measure section, including dropdowns and buttons.

        Returns
        -------
        html.Div
            A Dash HTML Div containing the Measure section components.
        """
        return html.Div([
            dcc.Store(id=self.plot_ready_store_id, data=False),  # Track if plot is ready
            html.Div([
                self.get_measure_dropdown(),
                self.get_xaxis_dropdown(),
                html.Button("Generate Plot", id=self.plot_button_id, n_clicks=0, style={'height': '36px'})
            ], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}),
            html.Div([
                dcc.Input(
                    id=self.filename_input_id,
                    type="text",
                    placeholder="Enter filename",
                    value="output.csv",
                    style={'margin-right': '10px', 'height': '36px', 'width': '200px'}
                ),
                html.Button("Save Data", id=self.save_button_id, n_clicks=0, disabled=True, style={'height': '36px', 'background-color': 'grey', 'color': 'white'})
            ], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center', 'margin-top': '20px'})
        ])

    def update_callbacks(self, callback_func, save_func):
        """
        Updates the callbacks for the Measure section components.

        Parameters
        ----------
        callback_func : callable
            A function to generate the plot based on the selected measure and x-axis.
        save_func : callable
            A function to save the data to a file when the "Save Data" button is clicked.
        """
        @self.app.callback(
            [Output("plot-image", "src"), Output(self.plot_ready_store_id, "data")],
            Input(self.plot_button_id, "n_clicks"),
            State(self.dropdown_id, "value"),
            State(self.xaxis_input_id, "value")
        )
        def trigger_callback(n_clicks, measure, xaxis):
            if n_clicks > 0:
                plot_src = callback_func(measure, xaxis)
                return plot_src, True
            return None, False

        @self.app.callback(
            [Output(self.save_button_id, "disabled"), Output(self.save_button_id, "style")],
            Input(self.plot_ready_store_id, "data")
        )
        def update_save_button_style(plot_ready):
            if plot_ready:
                return False, {'height': '36px', 'background-color': '#28a745', 'color': 'white'}  # Enabled (green button)
            return True, {'height': '36px', 'background-color': 'grey', 'color': 'white'}  # Disabled (grey button)

        @self.app.callback(
            Output(self.filename_input_id, "value"),
            Input(self.save_button_id, "n_clicks"),
            State(self.plot_ready_store_id, "data"),
            State(self.filename_input_id, "value"),
            State(self.dropdown_id, "value"),
            State(self.xaxis_input_id, "value")
        )
        def save_data(n_clicks, plot_ready, filename, measure, xaxis):
            if n_clicks > 0 and plot_ready:
                save_func(filename, measure, xaxis)
            return filename
