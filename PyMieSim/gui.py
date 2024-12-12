from dash import Dash, html, dcc, Input, Output, State
import matplotlib.pyplot as plt
import io
import base64
import webbrowser



class Section:
    def __init__(self, title, color, dropdown_options=None, dropdown_id=None, inputs=None, default_value=None):
        self.title = title
        self.color = color
        self.dropdown_options = dropdown_options
        self.dropdown_id = dropdown_id
        self.inputs = inputs
        self.default_value = default_value

    def create(self):
        """Create the section layout."""
        elements = [
            html.H2(self.title, style={'color': self.color}),
        ]

        if self.dropdown_options and self.dropdown_id:
            elements.append(html.Label(f"Select {self.title.lower()} type:"))
            elements.append(dcc.Dropdown(
                options=self.dropdown_options,
                value=self.default_value,  # Set default value
                placeholder=f"Choose a {self.title.lower()}...",
                id=self.dropdown_id
            ))
            elements.append(html.Br())

        if self.inputs:
            input_elements = []
            for input_id, input_placeholder in self.inputs:
                input_elements.append(html.Div([
                    html.Label(input_placeholder, style={'margin-right': '10px'}),
                    dcc.Input(
                        id=input_id,
                        type='text',
                        placeholder=input_placeholder,
                        style={'width': '150px'}
                    )
                ], style={'display': 'flex', 'align-items': 'center', 'margin-bottom': '10px'}))
            elements.extend(input_elements)

        return html.Div(elements, style={'border': '1px solid black', 'padding': '10px', 'margin': '10px'})


class OpticalSetupGUI:
    def __init__(self):
        """Initialize the Dash app and other settings."""
        self.app = Dash(__name__)
        self.setup_layout()
        self.setup_callbacks()

    def create_plot(self, data):
        """Generate a matplotlib plot and return it as a base64 image."""
        plt.switch_backend('Agg')  # Ensure Matplotlib doesn't attempt to open a GUI window

        plt.figure()
        plt.plot([0, 1, 2, 3], [0, 1, 4, 9], label=f"Plot for {data.get('source', {}).get('type', 'Unknown')}")
        plt.xlabel('X-axis')
        plt.ylabel('Y-axis')
        plt.title('Generated Plot')
        plt.legend()

        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        encoded_image = base64.b64encode(buf.read()).decode('utf-8')
        buf.close()
        return f"data:image/png;base64,{encoded_image}"

    def setup_layout(self):
        """Define the app layout."""
        source_section = Section(
            title="Source",
            color="blue",
            dropdown_options=[
                {'label': 'Laser', 'value': 'laser'},
                {'label': 'LED', 'value': 'led'},
                {'label': 'Sunlight', 'value': 'sunlight'}
            ],
            dropdown_id='source-dropdown',
            default_value='laser',
            inputs=[
                ('source-wavelength', "Wavelength (nm)"),
                ('source-power', "Power (W)"),
                ('source-na', "Numerical aperture (NA)")
            ]
        )

        scatterer_section = Section(
            title="Scatterer",
            color="green",
            dropdown_options=[
                {'label': 'Sphere', 'value': 'sphere'},
                {'label': 'Cylinder', 'value': 'cylinder'},
                {'label': 'CoreShell', 'value': 'coreshell'}
            ],
            dropdown_id='scatterer-dropdown',
            default_value='sphere',
            inputs=[
                ('scatterer-size', "Size (nm)"),
                ('scatterer-refractive-index', "Refractive index"),
                ('medium-refractive-index', "Medium refractive index")
            ]
        )

        detector_section = Section(
            title="Detector",
            color="red",
            dropdown_options=[
                {'label': 'Photodiode', 'value': 'photodiode'},
                {'label': 'Coherent', 'value': 'coherent'}
            ],
            dropdown_id='detector-dropdown',
            default_value='photodiode',
            inputs=[
                ('detector-na', "Numerical aperture"),
                ('detector-phi-offset', "Phi offset"),
                ('detector-gamma-offset', "Gamma offset"),
                ('detector-sampling', "Sampling"),
                ('detector-polarization_filter', "Polarization Filter")
            ]
        )

        self.app.layout = html.Div([
            html.Div([
                html.H1("Optical Simulation Interface", style={'text-align': 'center'}),
                source_section.create(),
                scatterer_section.create(),
                detector_section.create(),
                html.Button("Generate Plot", id="generate-plot", n_clicks=0)
            ], style={'width': '40%', 'float': 'left', 'padding': '10px'}),

            html.Div([
                html.H1("Visualization", style={'text-align': 'center'}),
                html.Img(id="plot-image", style={'width': '100%', 'height': 'auto'})
            ], style={'width': '50%', 'float': 'right', 'padding': '10px', 'border-left': '1px solid black'})
        ])

    def setup_callbacks(self):
        """Set up Dash callbacks."""
        @self.app.callback(
            Output("plot-image", "src"),
            Input("generate-plot", "n_clicks"),
            State("source-dropdown", "value"),
            State("source-wavelength", "value"),
            State("source-power", "value"),
            State("source-na", "value"),
            State("scatterer-dropdown", "value"),
            State("scatterer-size", "value"),
            State("scatterer-refractive-index", "value"),
            State("medium-refractive-index", "value"),
            State("detector-dropdown", "value"),
            State("detector-na", "value"),
            State("detector-phi-offset", "value"),
            State("detector-gamma-offset", "value"),
            State("detector-sampling", "value"),
            State("detector-polarization_filter", "value")
        )
        def update_plot(
            n_clicks, source_type, wavelength, power, na,
            scatter_type, size, refractive_index, medium_index,
            detector_type, detector_na, phi_offset, gamma_offset, sampling, polarization_filter
        ):
            if n_clicks > 0:
                data = {
                    'source': {
                        'type': source_type,
                        'wavelength': wavelength,
                        'power': power,
                        'numerical_aperture': na
                    },
                    'scatterer': {
                        'type': scatter_type,
                        'size': size,
                        'refractive_index': refractive_index,
                        'medium_refractive_index': medium_index
                    },
                    'detector': {
                        'type': detector_type,
                        'numerical_aperture': detector_na,
                        'phi_offset': phi_offset,
                        'gamma_offset': gamma_offset,
                        'sampling': sampling,
                        'polarization_filter': polarization_filter
                    }
                }
                return self.create_plot(data)
            return None

    def run(self):
        """Run the Dash app."""
        webbrowser.open("http://127.0.0.1:8050/", new=2)
        self.app.run_server(debug=True)



if __name__ == '__main__':
    gui = OpticalSetupGUI()
    gui.run()
