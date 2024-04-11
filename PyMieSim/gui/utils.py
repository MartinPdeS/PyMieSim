#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn, Dict
import numpy
import tkinter


class InputWidget:
    """
    A Widget class that encapsulates a GUI widget with specific properties.

    Attributes:
        default_value (float | str): The default value for the widget.
        label (str): A label for the widget used for identification.
        component_label (str): A label for the component part of the widget.
        multiplicative_factor (float | None): An optional factor by which the widget's value is multiplied.
        to_float (bool): A flag indicating whether the input should be converted to float. Defaults to True.
        is_permanent (bool): A flag indicating if the widget's value is permanent. Defaults to False.
        is_mappable (bool): A flag indicating if the widget's value can be mapped to other values. Defaults to True.

    Methods:
        update(): Updates the widget's value based on user input.
        get_input(): Retrieves the current input from the widget.
        process_input(): Processes the user input, converting it into a float or numpy array as appropriate.
    """

    def __init__(
            self,
            default_value: float | str,
            label: str,
            component_label: str,
            multiplicative_factor: float | None = None,
            to_float: bool = True,
            to_int: bool = False,
            is_permanent: bool = False,
            is_mappable: bool = True) -> None:
        """
        Initializes a new instance of the Widget class.
        """
        self.default_value = default_value
        self.tk_widget = tkinter.StringVar(value=str(default_value))
        self.label = label
        self.component_label = component_label
        self.to_float = to_float
        self.to_int = to_int
        self.value = None
        self.multiplicative_factor = multiplicative_factor
        self.is_permanent = is_permanent
        self.is_mappable = is_mappable
        self.update()

    def __repr__(self) -> str:
        return f"Widget(label={self.label})"

    def update(self) -> None:
        """
        Updates the widget's value based on the current user input.
        """
        self.value = self.process_input()

    def get_input(self) -> str:
        """
        Retrieves the current input from the tkinter StringVar associated with the widget.

        Returns:
            str: The current input value as a string.
        """
        return self.tk_widget.get()

    def process_input(self) -> numpy.ndarray | float:
        """
        Processes the user input, converting it into a float or numpy array based on the input format.

        Returns:
            numpy.ndarray | float: The processed input value, either as a float or numpy array.
        """
        user_input = self.get_input()
        values = numpy.nan  # Default case

        # Handling different input formats
        if "," in user_input:
            parts = [numpy.nan if p.lower() == 'none' else p.strip() for p in user_input.split(',')]
            values = numpy.asarray(parts)
        elif ":" in user_input:
            start, end, points = map(float, user_input.split(':'))
            values = numpy.linspace(start, end, int(points))
        else:
            values = numpy.nan if user_input.lower() == 'none' else user_input
            values = numpy.asarray(values)

        # Convert to float if necessary
        if self.to_float:
            values = values.astype(float)

        # Apply multiplicative factor if present
        if self.multiplicative_factor is not None:
            values *= self.multiplicative_factor

        if self.to_int:
            values = int(values)

        return values


class WidgetCollection:
    """
    A collection class for managing multiple Widget instances.

    This class facilitates the grouping of Widget instances, allowing for collective operations like
    updating widget values, clearing non-permanent widgets, and setting up widgets within a tkinter frame.

    Attributes:
        widgets (tuple[Widget, ...]): A tuple of Widget instances included in the collection.
    """

    def __init__(self, *widgets: InputWidget) -> None:
        """
        Initializes a new instance of WidgetCollection with a given set of Widget instances.

        Args:
            *widgets (Widget): Variable length Widget instances to be included in the collection.
        """
        self.widgets = widgets

    def to_component_dict(self) -> Dict[str, float | str]:
        """
        Creates a dictionary mapping component labels to their respective widget values.

        Returns:
            dict[str, float | str]: A dictionary where keys are component labels and values are widget values.
        """
        return {widget.component_label: widget.value for widget in self.widgets}

    def __getitem__(self, component_label: str) -> InputWidget:
        """
        Allows direct access to a Widget instance in the collection by its component label.

        Args:
            component_label (str): The component label of the desired Widget.

        Returns:
            Widget: The Widget instance with the matching component label, if found.
        """
        return next((widget for widget in self.widgets if widget.component_label == component_label), None)

    def clear_widgets(self) -> NoReturn:
        """
        Clears all non-permanent widgets from the tkinter frame.
        """
        for widget in self.widgets:
            if not widget.is_permanent:
                widget._label.destroy()
                widget._button.destroy()

    def update(self) -> NoReturn:
        """
        Updates the value of all widgets in the collection.
        """
        for widget in self.widgets:
            widget.update()

    def setup_widgets(self, frame: tkinter.Frame) -> NoReturn:
        """
        Sets up and packs the widgets within a specified tkinter frame.

        Args:
            frame (tkinter.Frame): The tkinter frame where widgets are to be packed.
        """
        for widget in self.widgets:
            widget._label = tkinter.Label(frame, text=widget.label)
            widget._label.pack(side=tkinter.BOTTOM)
            widget._button = tkinter.Entry(frame, textvariable=widget.tk_widget)
            widget._button.pack(side=tkinter.BOTTOM)
            # widget._button.configure(text='Password :')

    def __repr__(self) -> str:
        """
        Returns a string representation of the WidgetCollection instance, listing all widget labels.

        Returns:
            str: A space-separated string of widget labels.
        """
        return " ".join(str(widget) for widget in self.widgets)
