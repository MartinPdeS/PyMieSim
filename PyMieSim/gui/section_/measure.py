from dash import html, dcc, State, Input, Output
from PyMieSim.gui.helper import parse_string_to_array_or_float
import numpy

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
        self.download_id = "download-data"

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
            A Dash dropdown component populated dynamically with x-axis options.
        """
        return dcc.Dropdown(
            id=self.xaxis_input_id,
            options=[],
            value=None,
            style={'margin-right': '10px', 'height': '36px', 'width': '300px'}
        )

    def get_filename_input(self) -> dcc.Input:
        return dcc.Input(
            id=self.filename_input_id,
            type="text",
            placeholder="Enter filename",
            value="output.csv",
            style={'margin-right': '10px', 'height': '36px', 'width': '200px'}
        )

    def get_plot_button(self) -> html.Button:
        return html.Button("Generate Plot", id=self.plot_button_id, n_clicks=0, disabled=True, style={'height': '36px', 'background-color': 'grey', 'color': 'white'})

    def get_save_button(self) -> html.Button:
        return html.Button("Save Data", id=self.save_button_id, n_clicks=0, disabled=True, style={'height': '36px', 'background-color': 'grey', 'color': 'white'})

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
            dcc.Download(id=self.download_id),
            html.Div(
                [self.get_measure_dropdown(), self.get_xaxis_dropdown(), self.get_plot_button()],
                style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}
            ),
            html.Div(
                [self.get_filename_input(), self.get_save_button()],
                style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center', 'margin-top': '20px'}
            )
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
        def button_style(enabled: bool) -> tuple:
            """
            Helper function to generate button styles and states.

            Parameters
            ----------
            enabled : bool
                Whether the button should be enabled.

            Returns
            -------
            tuple
                A tuple containing (disabled, style) for the button.
            """
            if enabled:
                return False, {'height': '36px', 'background-color': '#28a745', 'color': 'white'}  # Green
            return True, {'height': '36px', 'background-color': 'grey', 'color': 'white'}  # Grey

        @self.app.callback(
            [Output("plot-image", "src"), Output(self.plot_ready_store_id, "data")],
            Input(self.plot_button_id, "n_clicks"),
            State(self.dropdown_id, "value"),
            State(self.xaxis_input_id, "value")
        )
        def trigger_callback(n_clicks: int, measure: str, xaxis: str) -> tuple:
            """
            Trigger the plot generation callback.

            Parameters
            ----------
            n_clicks : int
                Number of button clicks.
            measure : str
                Selected measure.
            xaxis : str
                Selected x-axis parameter.

            Returns
            -------
            tuple
                The plot source and a flag indicating plot readiness.
            """
            if n_clicks > 0:
                plot_src = callback_func(measure, xaxis)
                return plot_src, True
            return None, False

        @self.app.callback(
            [Output(self.save_button_id, "disabled"), Output(self.save_button_id, "style")],
            Input(self.plot_ready_store_id, "data")
        )
        def update_save_button_style(plot_ready: bool) -> tuple:
            """
            Update the save button's style based on plot readiness.

            Parameters
            ----------
            plot_ready : bool
                Whether the plot is ready.

            Returns
            -------
            tuple
                Disabled state and style for the save button.
            """
            return button_style(plot_ready)

        @self.app.callback(
            [Output(self.plot_button_id, "disabled"), Output(self.plot_button_id, "style")],
            Input(self.xaxis_input_id, "value")
        )
        def update_plot_button_style(xaxis_value: str) -> tuple:
            """
            Update the plot button's style based on the x-axis selection.

            Parameters
            ----------
            xaxis_value : str
                The selected x-axis value.

            Returns
            -------
            tuple
                Disabled state and style for the plot button.
            """
            return button_style(bool(xaxis_value))

        @self.app.callback(
            Output(self.download_id, "data"),
            Input(self.save_button_id, "n_clicks"),
            State(self.plot_ready_store_id, "data"),
            State(self.filename_input_id, "value"),
            State(self.dropdown_id, "value")
        )
        def save_data(n_clicks: int, plot_ready: bool, filename: str, measure: str) -> dict:
            """
            Trigger file download when the save button is clicked.

            Parameters
            ----------
            n_clicks : int
                Number of times the save button has been clicked.
            plot_ready : bool
                Whether the plot is ready.
            filename : str
                The file name for the downloaded data.
            measure : str
                Selected measure.
            xaxis : str
                Selected x-axis parameter.

            Returns
            -------
            dict or None
                File content and metadata for download, or None if conditions are not met.
            """
            if n_clicks > 0 and plot_ready:
                dataframe = save_func(filename=filename, measure=measure)
                content = dataframe.to_csv(index=False)
                return {"content": content, "filename": filename, "type": "text/csv"}
            return None


        @self.app.callback(
            Output(self.xaxis_input_id, "options"),
            Input(f"scatterer-sphere-data", "children"),
            Input(f"scatterer-coreshell-data", "children"),
            Input(f"scatterer-cylinder-data", "children"),
            Input(f"source-planewave-data", "children"),
            Input(f"source-gaussian-data", "children"),
            Input(f"detector-photodiode-data", "children"),
            Input(f"detector-coherentmode-data", "children"),
        )
        def update_xaxis_options(sphere_data, coreshell_data, cylinder_data, planewave_data, gaussian_data, photodiode_data, coherent_mode_data):
            """
            Dynamically update x-axis options based on Scatterer, Source, and Detector sections.

            Parameters
            ----------
            sphere_data : str
                The data from the sphere subsection of the scatterer.
            coreshell_data : str
                The data from the core-shell subsection of the scatterer.
            source_data : str
                The data from the source section.
            detector_data : str
                The data from the detector section.

            Returns
            -------
            list
                A list of dictionaries for the updated dropdown options.
            """
            options = []

            # Recompute _xaxis_options for all sections dynamically
            for section in [self.scatterer_section, self.source_section, self.detector_section]:
                section._xaxis_options = []  # Reset options
                section._xaxis_options_length = []  # Reset lengths

                # Parse inputs and update _xaxis_options
                for key, value in section.data.items():
                    try:
                        parsed_value = parse_string_to_array_or_float(value)
                        if isinstance(parsed_value, numpy.ndarray) and parsed_value.size > 1:
                            section._xaxis_options.append(f"{section.name}:{key}")
                            section._xaxis_options_length.append(parsed_value.size)
                    except ValueError:
                        pass  # Ignore invalid inputs

                # Add the recomputed options to the x-axis dropdown
                options.extend(
                    [{"label": f"{opt} ({size})", "value": opt} for opt, size in zip(section._xaxis_options, section._xaxis_options_length)]
                )

            return options
