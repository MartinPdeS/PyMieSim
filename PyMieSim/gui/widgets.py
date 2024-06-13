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

    def get_value(self):
        self.update()

        return self.value

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


class RadioButtonWidget(BaseWidget):

    def __init__(self, option_text: list, options_values: list, **kwargs):
        super().__init__(**kwargs)
        self.option_text = option_text
        self.tk_variable = tkinter.IntVar()
        self.options_values = options_values

    def update(self):
        pass

    def setup(self, row: int):
        self.tk_label = tkinter.Label(self.frame, text='Coupling mode: ')
        self.tk_label.grid(row=row, column=0, sticky="W")

        self.tk_widgets = []

        for column, text in enumerate(self.option_text):
            option = tkinter.Radiobutton(
                self.frame,
                text=text,
                variable=self.tk_variable,
                value=column,
            )

            option.grid(row=row, column=column + 1)
            self.tk_widgets.append(option)

    def destroy(self) -> NoReturn:
        for widget in self.tk_widgets:
            widget.destroy()

    def get_value(self):
        return self.options_values[self.tk_variable.get()]


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

    def get_value(self):
        self.update()

        return self.value

    def update(self) -> None:
        """
        Processes the user input, converting it into a float or numpy array based on the input format.

        Returns:
            numpy.ndarray | float: The processed input value, either as a float or numpy array.
        """
        user_input = self.tk_widget.get()
        value = numpy.nan  # Default case

        # Handling different input formats
        if "," in user_input:
            value = [numpy.nan if p.lower() == 'none' else p.strip() for p in user_input.split(',')]
        elif ":" in user_input:
            start, end, points = map(float, user_input.split(':'))
            value = numpy.linspace(start, end, int(points))
        else:
            value = numpy.nan if user_input.lower() == 'none' else user_input

        value = numpy.asarray(value)

        if self.dtype:
            value = value.astype(self.dtype)
        if self.multiplicative_factor is not None:
            value *= self.multiplicative_factor

        self.value = value



