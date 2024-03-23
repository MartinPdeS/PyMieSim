
from typing import NoReturn
from tkinter import ttk
import numpy as np
import tkinter


class Widget():
    def __init__(self, default_value, label: str, component_label: str, multiplicative_factor=None):
        self.default_value = default_value
        self.tk_widget = tkinter.StringVar(value=str(default_value))
        self.label = label
        self.component_label = component_label

    @property
    def user_input(self):
        return self.tk_widget.get()


class BaseTab:
    """
    A base class for tabs in the application's notebook widget, providing a structured way
    to create and populate tabs with GUI elements.

    Attributes:
        label (str): The label of the tab to be displayed.
        frame (ttk.Frame): The frame that holds the tab's content.
    """

    def __init__(self, notebook, label: str, main_window=None):
        """
        Initializes the tab and binds it to a notebook widget.

        Parameters:
            notebook (ttk.Notebook): The notebook widget to which this tab will be added.
            label (str): The label for the tab.
        """
        self.label = label
        self.frame = ttk.Frame(notebook)
        self.main_window = main_window
        notebook.add(self.frame, text=label)
        self.setup_tab()

    def parse_input(self, input_str: str, factor: float = None) -> float | np.ndarray:
        """
        Parses input strings for numerical values, supporting both single values and ranges.

        Parameters:
            input_str (str): The input string to parse.
            factor (float, optional): A factor to multiply with the parsed numbers, useful for unit conversion.

        Returns:
            float | np.ndarray: The parsed number(s) as a single float or an array of floats.
        """
        if ":" in input_str:
            parts = input_str.split(':')
            start, end, points = map(float, parts)
            values = np.linspace(start, end, int(points))

        elif "," in input_str:
            parts = input_str.split(',')
            values = [float(value) for value in parts]
            values = np.asarray(values)

        else:
            values = np.asarray(float(input_str))

        if factor is not None:
            values *= factor

        return values

    def setup_tab(self) -> NoReturn:
        """
        Sets up the tab's layout and widgets. To be implemented by subclasses.
        """
        raise NotImplementedError

    def get_inputs(self) -> list:
        user_inputs = []
        for key, val in self.variables.items():
            user_input = val['user_input']
            user_inputs.append(user_input.get())

        return user_inputs

    def update_user_input(self) -> list:
        for key, val in self.variables.items():
            user_input = val['user_input'].get()
            if val.get('to_float', True):
                user_input = self.parse_input(input_str=user_input, factor=val['factor'])

            self.variables[key]['values'] = user_input
