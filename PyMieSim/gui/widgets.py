#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Union, NoReturn, List
import numpy
import tkinter
from tkinter.ttk import Button

from PyMieSim.gui.singleton import datashelf

from pydantic.dataclasses import dataclass
from pydantic import ConfigDict


@dataclass(kw_only=True, config=ConfigDict(arbitrary_types_allowed=True))
class BaseWidget():
    label: str
    component_label: str
    dtype: type = None
    is_permanent: bool = False

    def __repr__(self) -> str:
        return f"Widget(label={self.label})"

    def destroy(self) -> NoReturn:
        self.tk_label.destroy()
        self.tk_widget.destroy()


@dataclass(kw_only=True, config=ConfigDict(arbitrary_types_allowed=True))
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
    default_options: int = 0
    options: List[str]
    value = None

    def initialize(self):
        """
        Empty function to allow uniform widget_collection.add_widgets accross all widgets
        """
        pass

    def setup(self, row: int = 0):
        # Sets up an attribute to datashelf that will hold the curent value of the combobox, for easy acces
        setattr(datashelf, f"{self.component_label}_selection", tkinter.StringVar(value=self.options[self.default_options]))
        self.text_variable = getattr(datashelf, f"{self.component_label}_selection")

        # Sets up the Combobox
        self.tk_label = tkinter.Label(self.frame, text=self.label)
        self.tk_widget = tkinter.ttk.Combobox(self.frame, values=self.options, textvariable=self.text_variable)
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


@dataclass(kw_only=True, config=ConfigDict(arbitrary_types_allowed=True))
class RadioButtonWidget(BaseWidget):

    option_text: list
    options_values: list
    can_be_axis: bool = False

    def initialize(self):
        self.tk_variable = tkinter.IntVar()

    def update(self):
        pass

    def setup(self, row: int):
        row = row + 1
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

            option.grid(row=row, sticky="W", column=column + 2)
            self.tk_widgets.append(option)

    def destroy(self) -> NoReturn:
        for widget in self.tk_widgets:
            widget.destroy()
        self.tk_label.destroy()

    def get_value(self):
        return self.options_values[self.tk_variable.get()]


@dataclass(kw_only=True, config=ConfigDict(arbitrary_types_allowed=True))
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
    default_value: Union[float | str]
    multiplicative_factor: float | None = None
    can_be_axis: bool = True

    def initialize(self) -> NoReturn:
        """
        Initializes a new instance of the Widget class.
        """
        self.lol = tkinter.StringVar(value=str(self.default_value))
        self.update()

    def setup(self, row: int):
        row += 1
        self.tk_label = tkinter.Label(self.frame, text=self.label)
        self.tk_label.grid(row=row + 1, column=0, sticky="W", pady=2)
        self.tk_widget = tkinter.Entry(self.frame, textvariable=self.lol)
        self.tk_widget.grid(row=row + 1, column=1, sticky="W", pady=2)

        # Adds the radiobuttons used to select wether this variable is used as an axis
        if self.can_be_axis:
            self.tk_radio_button_1 = tkinter.Radiobutton(self.frame, variable=datashelf.x_axis_label_widget, value=self.component_label)
            self.tk_radio_button_1.grid(row=row + 1, column=2, sticky="W", pady=2)
            self.tk_radio_button_2 = tkinter.Radiobutton(self.frame, variable=datashelf.STD_axis_label_widget, value=self.component_label)
            self.tk_radio_button_2.grid(row=row + 1, column=3, sticky="W", pady=2)

    def get_value(self):
        self.update()

        return self.value

    def update(self) -> None:
        """
        Processes the user input, converting it into a float or numpy array based on the input format.

        Returns:
            numpy.ndarray | float: The processed input value, either as a float or numpy array.
        """
        user_input = self.lol.get()
        value = numpy.nan  # Default case

        # Handling different input formats
        if "," in user_input:
            value = [numpy.nan if p.lower() == 'none' else p.strip() for p in user_input.split(',')]
        elif ":" in user_input:
            start, end, points = map(float, user_input.split(':'))
            value = numpy.linspace(start, end, int(points))
        else:
            value = numpy.nan if user_input.lower() == 'none' else user_input
            value = numpy.atleast_1d(value)

        value = numpy.asarray(value)

        if self.dtype:
            value = value.astype(self.dtype)

        if self.multiplicative_factor is not None:
            value *= self.multiplicative_factor

        self.value = value

    def destroy(self) -> NoReturn:
        self.tk_label.destroy()
        self.tk_widget.destroy()
        if self.can_be_axis:
            self.tk_radio_button_1.destroy()
            self.tk_radio_button_2.destroy()


class ControlWidget():
    """
    Buttons that will allow the user to calculate the simulation, save it,
    export it and reset their std-axis selection.

    Attribute:
    frame (ttk.Frame): given in WidgetCollection
    command (method from ControlTab): given in WidgetCollection
    component_label (str): given in WidgetCollection
    label (str): the name of the button
    style (str): the style of the button
    """

    def __init__(self, label: str, style: str = "Large.TButton") -> None:
        self.label = label
        self.style = style

    def __repr__(self) -> str:
        return f"Widget(label={self.label})"

    def setup(self, column):
        self.button = Button(
            self.frame,
            text=self.label,
            style=self.style,
            command=self.command
        )
        self.button.grid(row=0, column=column, sticky='ew')

    def invoke(self):
        self.button.invoke()

# -
