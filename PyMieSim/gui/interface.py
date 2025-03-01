from dash import Dash, html, dcc
import matplotlib.pyplot as plt
import io
import base64
import webbrowser
from PyMieSim.gui.section import SourceSection, ScattererSection, DetectorSection, MeasureSection
from PyMieSim.gui.helper import interface

dcc_store_id = "input-store"

class OpticalSetupGUI:
    """
    A class to manage the overall GUI for the optical simulation interface.

    This class integrates various sections, including Source, Scatterer, Detector, and Measure sections,
    and sets up the layout and callbacks for the Dash application.
    """

    def __init__(self):
        """
        Initialize the Dash app and other settings.

        This method sets up the Dash application, initializes the individual sections, and prepares the layout
        and callbacks.
        """
        plt.switch_backend('Agg')
        self.app = Dash(__name__, external_stylesheets=["https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css"])

        self.source_section = SourceSection(self.app)
        self.scatterer_section = ScattererSection(self.app)
        self.detector_section = DetectorSection(self.app)
        self.measure_section = MeasureSection(self.app, self.scatterer_section, self.source_section, self.detector_section)

        self.setup_layout()
        self.setup_callbacks()

    def create_plot(self, measure: str, xaxis: str):
        """
        Generate a matplotlib plot and return it as a base64 image.

        Parameters
        ----------
        measure : str
            The type of measure to plot (e.g., Qsca, Qext).
        xaxis : str
            The parameter to plot on the x-axis.

        Returns
        -------
        str
            A base64-encoded string of the generated plot image.
        """
        data = interface(
            source_kwargs=self.source_section.data,
            scatterer_kwargs=self.scatterer_section.data,
            detector_kwargs=self.detector_section.data,
            measure=measure
        )

        data.plot(x=xaxis)

        buf = io.BytesIO()
        plt.savefig(buf, format='png')
        buf.seek(0)
        encoded_image = base64.b64encode(buf.read()).decode('utf-8')
        buf.close()
        return f"data:image/png;base64,{encoded_image}"

    def save_func(self, filename: str, measure: str):
        """
        Save the simulation data to a CSV file.

        Parameters
        ----------
        filename : str
            The name of the file to save the data.
        measure : str
            The type of measure to save.
        xaxis : str
            The parameter to include on the x-axis in the saved data.
        """
        return interface(
            source_kwargs=self.source_section.data,
            scatterer_kwargs=self.scatterer_section.data,
            detector_kwargs=self.detector_section.data,
            measure=measure,
            add_units=False
        )

    def setup_layout(self):
        """
        Define the layout of the Dash app.

        This method organizes the sections and visualization components in the application layout.
        """
        self.app.layout = html.Div([
            dcc.Store(id=dcc_store_id),  # Store to save input values
            html.Div([
                html.H1("Optical Simulation Interface", style={'text-align': 'center'}),
                *self.source_section.create(),
                *self.scatterer_section.create(),
                *self.detector_section.create(),
            ], style={'width': '40%', 'float': 'left', 'padding': '10px'}),

            html.Div([
                html.H1("Visualization", style={'text-align': 'center'}),
                self.measure_section.create(),
                html.Img(id="plot-image", style={'width': '100%', 'height': 'auto', 'margin-top': '20px'})
            ], style={'width': '55%', 'float': 'right', 'padding': '10px', 'border-left': '1px solid black'})
        ])

    def setup_callbacks(self):
        """
        Set up Dash callbacks for the GUI.

        This method links user interactions with the relevant functions to update the plot and save data.
        """
        self.measure_section.update_callbacks(self.create_plot, self.save_func)

    def run(self, host: str = "0.0.0.0", port: str = "8050", open_browser: bool = False, debug: bool = False):
        """
        Run the Dash app.

        This method starts the Dash server and opens the application in the default web browser.
        """
        if open_browser:
            webbrowser.open(f"http://{host}:{port}/", new=2)
            self.app.run_server(debug=debug)

        else:
            self.app.run_server(debug=debug, host=host, port=port)


if __name__ == '__main__':
    _gui = OpticalSetupGUI()
    _gui.run(host='0.0.0.0', port='8050', open_browser=False, debug=True)