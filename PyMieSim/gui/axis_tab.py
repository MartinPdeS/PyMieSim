#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import tkinter
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.utils import Widget, WidgetCollection


class AxisTab(BaseTab):
    measure_map = scatterer.Sphere.available_measure_list

    def __init__(self, *args, other_tabs: list, **kwargs):
        """
        Initializes the Axis Configuration tab with predefined variables for selecting
        the axes for plotting the simulation results.
        """
        self.other_tabs = other_tabs

        super().__init__(*args, **kwargs)

    def clear_button(self) -> NoReturn:
        for element in self.non_permanent_widget:
            element.destroy()

    def setup(self) -> NoReturn:
        """
        Sets up the GUI elements for the Axis Configuration tab, including labels and comboboxes
        for selecting the x and y axis variables. This method provides a user interface for
        specifying the variables to be plotted on the axes.
        """
        self.x_axis_options = list(self.axis_mapping.keys())
        self.y_axis_options = list(self.measure_map.keys())

        self.widget_collection = WidgetCollection(
            Widget(default_value='phi', label='x-axis', component_label='x_axis', to_float=False),
            Widget(default_value='coupling', label='y-axis', component_label='y_axis', to_float=False),
            Widget(default_value='none', label='STD-axis', component_label='std_axis', to_float=False),
        )
        self.non_permanent_widget = []

        self.setup_combox(sub_text='STD axis', component_label='std_axis')
        self.setup_combox(sub_text='x axis', component_label='x_axis')
        self.setup_combox(sub_text='y axis', component_label='y_axis')

    def setup_combox(self, sub_text: str, component_label: str) -> NoReturn:
        label = tkinter.Label(
            self.frame,
            text=sub_text
        )

        label.pack(side=tkinter.BOTTOM)

        self.non_permanent_widget.append(label)

        combox = tkinter.ttk.Combobox(
            self.frame,
            textvariable=self.widget_collection[component_label].tk_widget,
            values=self.x_axis_options,
            state="readonly"
        )

        combox.pack(side=tkinter.BOTTOM)

        self.non_permanent_widget.append(combox)

    @property
    def x_axis(self) -> str:
        x_axis = self.widget_collection['x_axis'].tk_widget.get()

        return self.axis_mapping[x_axis]

    @property
    def std_axis(self) -> str:
        std_axis = self.widget_collection['std_axis'].tk_widget.get()

        return self.axis_mapping[std_axis]

    @property
    def axis_mapping(self) -> dict:
        _axis_mapping = {}
        for tab in self.other_tabs:
            _axis_mapping.update(tab.mapping)

        return _axis_mapping
