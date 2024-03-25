#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import numpy as np
import tkinter


def parse_input(input_str: str, factor: float = None) -> float | np.ndarray:
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


class Widget():
    def __init__(
            self,
            default_value,
            label: str,
            component_label: str,
            multiplicative_factor=None,
            to_float: bool = True,
            is_permanent: bool = False,
            is_mappable: bool = True):

        self.default_value = default_value
        self.tk_widget = tkinter.StringVar(value=str(default_value))
        self.label = label
        self.component_label = component_label
        self.to_float = to_float
        self.value = None
        self.multiplicative_factor = multiplicative_factor
        self.is_permanent = is_permanent
        self.is_mappable = is_mappable
        self.update()

    def __repr__(self):
        return self.label

    def update(self):
        self.user_input = _user_input = self.tk_widget.get()
        if self.to_float:
            _user_input = parse_input(input_str=_user_input, factor=self.multiplicative_factor)

        self.value = _user_input

    def get_input(self):
        return self.tk_widget.get()


class WidgetCollection():
    def __init__(self, *widgets):
        self.widgets = widgets

    def to_component_dict(self) -> dict:
        return {
            widget.component_label: widget.value for widget in self.widgets
        }

    def __getitem__(self, component_label: str) -> Widget:
        for widget in self.widgets:
            if widget.component_label == component_label:
                return widget

    def clear_widgets(self) -> NoReturn:
        for widget in self.widgets:
            if not widget.is_permanent:
                widget._label.destroy()
                widget._button.destroy()

    def update(self) -> NoReturn:
        for widget in self.widgets:
            widget.update()

    def setup_widgets(self, frame) -> NoReturn:
        for widget in self.widgets:
            widget._label = tkinter.Label(frame, text=widget.label)
            widget._label.pack(side=tkinter.BOTTOM)

            widget._button = tkinter.Entry(frame, textvariable=widget.tk_widget)
            widget._button.pack(side=tkinter.BOTTOM)

    def __repr__(self) -> str:
        return " ".join(str(widget) for widget in self.widgets)
