from dash import html, dcc, State, Input, Output
from PyMieSim.units import nanometer, degree, milliwatt


length_units = nanometer
power_units = milliwatt
angle_units = degree

class Section:
    def create_call_backs(self):
        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [Input(id, "value") for id in self.input_id]
        self.update_callbacks()

    def create(self):
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
            html.Div(id=f"{self.name}-dropdown-data", style={"display": "none"}) # Hidden placeholder for scatterer section
        )


class SourceSection(Section):
    name = 'source'
    def __init__(self, app):
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
        """Set up callbacks to update the data dictionary."""
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
    name = 'scatterer'
    def __init__(self, app):
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
                "default": "1.5"
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
        """Set up callbacks to update the data dictionary."""
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
    name = 'detector'
    def __init__(self, app):
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
        """Set up callbacks to update the data dictionary."""
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
    def __init__(self, app, scatterer_section):
        self.app = app
        self.scatterer_section = scatterer_section
        self.dropdown_id = "measure-input"
        self.button_id = "generate-plot"
        self.xaxis_input_id = "xaxis-input"
        self.data = "Qsca"  # Default measure value

    def get_measure_dropdown(self) -> dcc.Dropdown:
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
        return dcc.Dropdown(
            id=self.xaxis_input_id,
            options=[
                {"label": "scatterer:" + k, "value": "scatterer:" + k} for k in self.scatterer_section.inputs.keys()],
            value='1',  # Default value
            style={'margin-right': '10px', 'height': '36px', 'width': '300px'}
        )

    def create(self) -> html.Div:
        return html.Div([
            html.Div([
                self.get_measure_dropdown(), self.get_xaxis_dropdown(), html.Button("Generate Plot", id=self.button_id, n_clicks=0, style={'height': '36px'})],
                style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})
        ])

    def update_callbacks(self, callback_func):
        @self.app.callback(
            Output("plot-image", "src"),
            Input(self.button_id, "n_clicks"),
            State(self.dropdown_id, "value"),
            State(self.xaxis_input_id, "value")
        )
        def trigger_callback(n_clicks, measure, xaxis):
            if n_clicks > 0:
                return callback_func(measure, xaxis)
            return None
