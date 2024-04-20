#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import numpy
import tkinter


class BaseWidget():
    def __init__(self, label: str, component_label: str, dtype: type = None, is_permanent: bool = False) -> None:
        """
        Initializes a new instance of the Widget class.
        """
        self.label = label
        self.component_label = component_label
        self.dtype = dtype
        self.is_permanent = is_permanent

    def __repr__(self) -> str:
        return f"Widget(label={self.label})"

    def destroy(self) -> NoReturn:
        self.tk_label.destroy()
        self.tk_widget.destroy()


class ComBoxWidget(BaseWidget):
    """
    A Widget class that encapsulates a GUI widget with specific properties.

    Attributes:
        default_value (float | str): The default value for the widget.
        label (str): A label for the widget used for identification.
        component_label (str): A label for the component part of the widget.
        multiplicative_factor (float | None): An optional factor by which the widget's value is multiplied.
        to_float (bool): A flag indicating whether the input should be converted to float. Defaults to True.
        is_permanent (bool): A flag indicating if the widget's value is permanent. Defaults to False.

    Methods:
        update(): Updates the widget's value based on user input.
        get_input(): Retrieves the current input from the widget.
        process_input(): Processes the user input, converting it into a float or numpy array as appropriate.
    """

    def __init__(self, default_options: int = 0, options: list = [], **kwargs) -> None:
        """
        Initializes a new instance of the Widget class.
        """
        super().__init__(**kwargs)
        self.default_options = default_options
        self.value = None
        self.options = options

    def setup(self, row: int = 0):
        self.tk_label = tkinter.Label(self.frame, text=self.label)
        self.tk_widget = tkinter.ttk.Combobox(self.frame, values=self.options)
        self.tk_widget.current(self.default_options)

        self.tk_widget.grid(row=row, column=1)
        self.tk_label.grid(row=row, column=0)

    def update(self) -> None:
        """
        Updates the widget's value based on the current user input.
        """
        self.value = self.tk_widget.get()

    def get_input(self) -> str:
        """
        Retrieves the current input from the tkinter StringVar associated with the widget.

        Returns:
            str: The current input value as a string.
        """
        return self.tk_widget.get()


class InputWidget(BaseWidget):
    """
    A Widget class that encapsulates a GUI widget with specific properties.

    Attributes:
        default_value (float | str): The default value for the widget.
        label (str): A label for the widget used for identification.
        component_label (str): A label for the component part of the widget.
        multiplicative_factor (float | None): An optional factor by which the widget's value is multiplied.
        to_float (bool): A flag indicating whether the input should be converted to float. Defaults to True.
        is_permanent (bool): A flag indicating if the widget's value is permanent. Defaults to False.

    Methods:
        update(): Updates the widget's value based on user input.
        process_input(): Processes the user input, converting it into a float or numpy array as appropriate.
    """

    def __init__(self, default_value: float | str, multiplicative_factor: float | None = None, **kwargs) -> None:
        """
        Initializes a new instance of the Widget class.
        """
        super().__init__(**kwargs)
        self.default_value = default_value
        self.tk_widget = tkinter.StringVar(value=str(default_value))
        self.value = None
        self.multiplicative_factor = multiplicative_factor
        self.update()

    def setup(self, row: int):
        self.tk_label = tkinter.Label(self.frame, text=self.label)
        self.tk_label.grid(row=row + 1, column=0, sticky="W", pady=2)
        self.tk_widget = tkinter.Entry(self.frame, textvariable=self.tk_widget)
        self.tk_widget.grid(row=row + 1, column=1, sticky="W", pady=2)

    def update(self) -> None:
        """
        Processes the user input, converting it into a float or numpy array based on the input format.

        Returns:
            numpy.ndarray | float: The processed input value, either as a float or numpy array.
        """
        user_input = self.tk_widget.get()
        values = numpy.nan  # Default case

        # Handling different input formats
        if "," in user_input:
            values = self.process_coma_input(user_input=user_input)
        elif ":" in user_input:
            values = self.process_coma_double_point(user_input=user_input)
        else:
            values = self.process_coma_else(user_input=user_input)

        if self.dtype:
            values = values.astype(self.dtype)
        if self.multiplicative_factor is not None:
            values *= self.multiplicative_factor

        self.value = values

    def process_coma_input(self, user_input: str) -> numpy.ndarray | float:
        parts = [numpy.nan if p.lower() == 'none' else p.strip() for p in user_input.split(',')]
        return numpy.asarray(parts)

    def process_coma_double_point(self, user_input: str) -> numpy.ndarray | float:
        start, end, points = map(float, user_input.split(':'))
        return numpy.linspace(start, end, int(points))

    def process_coma_else(self, user_input: str) -> numpy.ndarray | float:
        values = numpy.nan if user_input.lower() == 'none' else user_input
        return numpy.asarray(values)
