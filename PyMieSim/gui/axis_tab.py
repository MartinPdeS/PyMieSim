#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.utils import WidgetCollection, ComBoxWidget


class AxisTab(BaseTab):
    measure_map = scatterer.Sphere.available_measure_list

    def __init__(self, *args, other_tabs: list, **kwargs):
        """
        Initializes the Axis Configuration tab with predefined variables for selecting
        the axes for plotting the simulation results.
        """
        self.other_tabs = other_tabs

        super().__init__(*args, **kwargs)

    def setup(self) -> NoReturn:
        """
        Sets up the GUI elements for the Axis Configuration tab, including labels and comboboxes
        for selecting the x and y axis variables. This method provides a user interface for
        specifying the variables to be plotted on the axes.
        """
        self.x_axis_options = list(self.axis_mapping.keys())
        self.y_axis_options = list(self.measure_map.keys())

        self.widget_collection = WidgetCollection(
            ComBoxWidget(label='x-axis', component_label='x_axis', options=self.x_axis_options, frame=self.frame, default_options=0),
            ComBoxWidget(label='y-axis', component_label='y_axis', options=self.y_axis_options, frame=self.frame, default_options=0),
            ComBoxWidget(label='STD-axis', component_label='std_axis', options=['none', *self.x_axis_options], frame=self.frame, default_options=0),
        )

        self.widget_collection.setup_widgets(frame=self.frame)

    @property
    def x_axis(self) -> str:
        x_axis = self.widget_collection.widgets[0].tk_widget.get()
        return self.axis_mapping[x_axis]

    @property
    def std_axis(self) -> str:
        std_axis = self.widget_collection.widgets[2].tk_widget.get()
        if std_axis == 'none':
            return None
        return self.axis_mapping[std_axis]

    @property
    def axis_mapping(self) -> dict:
        _axis_mapping = {}
        for tab in self.other_tabs:
            _axis_mapping.update(tab.mapping)

        return _axis_mapping
