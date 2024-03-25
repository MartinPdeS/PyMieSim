#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import numpy
import tkinter


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
        _user_input = self.process_input()

        self.value = _user_input

    def get_input(self):
        return self.tk_widget.get()

    def process_input(self):
        user_input = self.get_input()

        if "," in user_input:
            parts = user_input.split(',')
            values = [value.strip() for value in parts]
            values = numpy.asarray(values)

        elif ":" in user_input:
            start, end, points = user_input.split(':')
            values = numpy.linspace(float(start), float(end), int(points))

        else:
            values = numpy.asarray(user_input)

        if self.to_float:
            values = values.astype(float)

        if self.multiplicative_factor:
            values *= self.multiplicative_factor

        return values


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
