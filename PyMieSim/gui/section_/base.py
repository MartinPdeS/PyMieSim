from dash import html, dcc, Input, Output
from PyMieSim import units
from PyMieSim.gui.helper import parse_string_to_array_or_float
import numpy

length_units = units.nanometer
power_units = units.milliwatt
angle_units = units.degree

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
        title = html.H2(self.title, style={'color': self.color})
        type_dropdown = html.Div([
                html.Label(f"Select {self.title} Type:"),
                html.Div([dcc.Dropdown(id=self.dropdown_id, options=self.dropdown_options, value=self.default_value)])
            ],
            style={'margin-bottom': '20px'}
        )

        layout = [
            title,
            type_dropdown,
        ]

        for input in self.inputs.values():
            user_option = html.Div([
                html.Label(input['label'], style={'margin-right': '10px'}),
                dcc.Input(id=input['id'], type='text', value=input['default'], style={'width': '200px'})
                ],
                style={'margin-bottom': '10px'}
            )
            layout.append(user_option)

        return (
            html.Div(layout, style={'padding': '10px', 'border': '1px solid black', 'margin': '10px'}),
            html.Div(id=f"{self.name}-dropdown-data", style={"display": "none"})  # Hidden placeholder for data
        )

    def update_callbacks(self):
        """
        Set up callbacks to update the data dictionary and x-axis options.
        """
        @self.app.callback(Output(f"{self.dropdown_id}-data", "children"), *self.call_backs)
        def update_data(object_type, *inputs):
            """
            Update the data dictionary and x-axis options based on user inputs.

            Parameters
            ----------
            object_type : str
                The type of the source selected in the dropdown.
            inputs : list
                The values of the inputs from the source section.

            Returns
            -------
            str
                A message indicating the data has been updated.
            """
            # Map the input values to their corresponding keys
            input_values = dict(zip(self.inputs.keys(), inputs))

            # Parse inputs and update _xaxis_options if input is an array
            self._xaxis_options = []
            for key, value in input_values.items():
                try:
                    parsed_value = parse_string_to_array_or_float(value)
                    if isinstance(parsed_value, numpy.ndarray) and parsed_value.size > 1:
                        self._xaxis_options.append(f"{self.name}:{key}")
                        self._xaxis_options_length.append(len(parsed_value))
                except ValueError:
                    pass  # Ignore invalid inputs

            # Update the data dictionary
            self.data = {'type': object_type, **input_values}

            return "Data Updated"
