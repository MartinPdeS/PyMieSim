from dash import html, dcc, State, Input, Output
from PyMieSim.experiment.source.gaussian import Gaussian
from PyMieSim.experiment.scatterer.sphere import Sphere
from PyMieSim.experiment.detector.photodiode import Photodiode
from PyMieSim.experiment import Setup
from PyMieSim.units import nanometer, degree, milliwatt, AU, RIU
import numpy


length_units = nanometer
power_units = milliwatt
angle_units = degree


def parse_string_to_array_or_float(input_str):
    """
    Parse a string to return either a numpy array or a float.

    If the string is in the format 'start:end:count', return np.linspace(start, end, count).
    If the string is a single numeric value, return it as a float.
    """
    try:
        if ":" in input_str:
            start, end, count = map(float, input_str.split(":"))
            return numpy.linspace(start, end, int(count))
        else:
            return float(input_str)
    except ValueError:
        raise ValueError("Invalid input string format. Expected 'start:end:count' or a single numeric value.")

def interface(source_kwargs: dict, scatterer_kwargs: dict, detector_kwargs: dict, measure: str):
    source = Gaussian(
        wavelength=parse_string_to_array_or_float(source_kwargs['wavelength']) * length_units,
        polarization=parse_string_to_array_or_float(source_kwargs['polarization']) * angle_units,
        NA=parse_string_to_array_or_float(source_kwargs['NA']) * AU,
        optical_power=parse_string_to_array_or_float(source_kwargs['optical_power']) * milliwatt
    )

    scatterer = Sphere(
        diameter=parse_string_to_array_or_float(scatterer_kwargs['diameter']) * length_units,
        property=parse_string_to_array_or_float(scatterer_kwargs['property']) * RIU,
        medium_property=parse_string_to_array_or_float(scatterer_kwargs['medium_property']) * RIU,
        source=source
    )

    detector = Photodiode(
        NA=parse_string_to_array_or_float(detector_kwargs['NA']) * AU,
        gamma_offset=parse_string_to_array_or_float(detector_kwargs['gamma_offset']) * angle_units,
        phi_offset=parse_string_to_array_or_float(detector_kwargs['phi_offset']) * angle_units,
        sampling=parse_string_to_array_or_float(detector_kwargs['sampling']) * AU,
    )

    setup = Setup(source=source, scatterer=scatterer, detector=detector)

    dataframe = setup.get(measure)

    return dataframe

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

        for input_id, placeholder, default_value in zip(self.input_id, self.input_repr, self.inputs_default):
            layout.append(html.Div([
                html.Label(placeholder, style={'margin-right': '10px'}),
                dcc.Input(
                    id=input_id,
                    type='text',
                    value=default_value,
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

        self.input_id = (f'{self.name}-wavelength', f'{self.name}-optical_power', f'{self.name}-NA', f'{self.name}-polarization')
        self.input_repr = (f"Wavelength [{length_units}]", f"Power [{power_units}]", "Numerical aperture [NA]", f"Polarization [{angle_units}]")
        self.inputs_default = ('650', '5', '0.2', '0')

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [Input(id, "value") for id in self.input_id]

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

        self.input_id = (f'{self.name}-diameter', f'{self.name}-property', f'{self.name}-medium-property')
        self.input_repr = (f"Diameter [{length_units}]", "Scatterer Property", "Medium Property")
        self.inputs_default = ('100:20000:200', '1.5', '1.33')

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [Input(id, "value") for id in self.input_id]
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
        self.input_id = (f'{self.name}-na', f'{self.name}-phi-offset', f'{self.name}-gamma-offset', f'{self.name}-sampling', f'{self.name}-polarization_filter')

        self.input_repr = ("Numerical aperture", f"Phi offset [{angle_units}]", f"Gamma offset [{angle_units}]", "Sampling", f"Polarization Filter [{angle_units}]")

        self.inputs_default = ('0.9', '0', '0', '1000', 'None')

        self.call_backs = [Input(f"{self.name}-dropdown", "value")] + [Input(id, "value") for id in self.input_id]

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

    def create(self):
        return html.Div([
            html.Div([
                dcc.Dropdown(
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
                        {'label': 'g', 'value': 'g'}
                    ],
                    value='Qsca',
                    style={'margin-right': '10px', 'margin-bottom': '0px', 'height': '36px'}
                ),
                dcc.Input(
                    id=self.xaxis_input_id,
                    type="text",
                    placeholder="Enter X-axis",
                    value="scatterer:diameter",
                    style={'margin-right': '10px', 'height': '36px'}
                ),
                html.Button("Generate Plot", id=self.button_id, n_clicks=0, style={'height': '36px'})
            ], style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'})
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