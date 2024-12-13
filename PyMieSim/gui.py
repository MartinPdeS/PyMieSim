from dash import Dash, html, Input, Output, State, dcc
import matplotlib.pyplot as plt
import io
import base64
import webbrowser
from PyMieSim.gui_section import SourceSection, ScattererSection, DetectorSection, MeasureSection, interface

dcc_store_id = "input-store"

class OpticalSetupGUI:
    def __init__(self):
        """Initialize the Dash app and other settings."""
        plt.switch_backend('Agg')
        self.app = Dash(__name__, external_stylesheets=["https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css"])

        self.source_section = SourceSection(self.app)
        self.scatterer_section = ScattererSection(self.app)
        self.detector_section = DetectorSection(self.app)
        self.measure_section = MeasureSection(self.app, self.scatterer_section)

        self.setup_layout()
        self.setup_callbacks()

    def create_plot(self, measure: str, xaxis: str):
        """Generate a matplotlib plot and return it as a base64 image."""
        data = interface(
            source_kwargs=self.source_section.data,
            scatterer_kwargs=self.scatterer_section.data,
            detector_kwargs=self.detector_section.data,
            measure=measure
        )

        data.plot_data(x=xaxis)

        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        encoded_image = base64.b64encode(buf.read()).decode('utf-8')
        buf.close()
        return f"data:image/png;base64,{encoded_image}"

    def setup_layout(self):
        """Define the app layout."""
        self.app.layout = html.Div([
            dcc.Store(id=dcc_store_id),  # Store to save input values
            html.Div([
                html.H1("Optical Simulation Interface", style={'text-align': 'center'}),
                *self.source_section.create(),
                *self.scatterer_section.create(),
                *self.detector_section.create(),
            ], style={'width': '45%', 'float': 'left', 'padding': '10px'}),

            html.Div([
                html.H1("Visualization", style={'text-align': 'center'}),
                self.measure_section.create(),
                html.Img(id="plot-image", style={'width': '100%', 'height': 'auto', 'margin-top': '20px'})
            ], style={'width': '45%', 'float': 'right', 'padding': '10px', 'border-left': '1px solid black'})
        ])

    def setup_callbacks(self):
        """Set up Dash callbacks."""
        self.measure_section.update_callbacks(self.create_plot)

    def run(self):
        """Run the Dash app."""
        webbrowser.open("http://127.0.0.1:8050/", new=2)
        self.app.run_server(debug=True)

if __name__ == '__main__':
    gui = OpticalSetupGUI()
    gui.run()
