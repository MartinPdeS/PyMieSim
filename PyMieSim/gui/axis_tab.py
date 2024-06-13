#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn, Optional, Dict
from tkinter import ttk
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widgets import ComBoxWidget
from PyMieSim.gui.widget_collection import WidgetCollection


class AxisTab(BaseTab):
    measure_map = scatterer.Sphere.available_measure_list

    def __init__(self, master: ttk.Notebook, label: str, other_tabs: list[BaseTab], **kwargs) -> NoReturn:
        """
        Initializes the Axis Configuration tab with references to other tabs to gather possible axis choices.

        Args:
            master (ttk.Notebook): The notebook widget this tab will be part of.
            label (str): The label for the tab.
            other_tabs (List[BaseTab]): List of other tab instances to reference for setting up axis mappings.
        """
        self.other_tabs = other_tabs
        super().__init__(master, label=label)
        self.setup()

    def setup(self) -> NoReturn:
        """
        Sets up the UI elements for axis configuration, including comboboxes for selecting
        variables for the x-axis, y-axis, and an optional standard deviation (STD) axis.
        """
        self.x_axis_options = list(self.axis_mapping.keys())
        self.y_axis_options = list(self.measure_map.keys())

        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(
            ComBoxWidget(label='x-axis', component_label='x_axis', options=self.x_axis_options, default_options=11),
            ComBoxWidget(label='y-axis', component_label='y_axis', options=self.y_axis_options, default_options=21),
            ComBoxWidget(label='STD-axis', component_label='std_axis', options=['none', *self.x_axis_options], default_options=0),
        )

        self.widget_collection.setup_widgets()

    @property
    def x_axis(self) -> str:
        """
        Retrieves the selected x-axis variable from the widget collection.

        Returns:
            str: The key in the axis mapping corresponding to the selected x-axis variable.
        """
        x_axis = self.widget_collection.widgets[0].tk_widget.get()
        return self.axis_mapping[x_axis]

    @property
    def std_axis(self) -> Optional[str]:
        """
        Retrieves the selected standard deviation axis variable, if any.

        Returns:
            Optional[str]: The key in the axis mapping corresponding to the selected standard deviation axis variable,
            or None if 'none' is selected.
        """
        std_axis = self.widget_collection.widgets[2].tk_widget.get()
        if std_axis == 'none':
            return None
        return self.axis_mapping[std_axis]

    @property
    def axis_mapping(self) -> Dict[str, str]:
        """
        Combines mappings from all other tabs to provide a comprehensive dictionary of available axis options.

        Returns:
            Dict[str, str]: A dictionary mapping UI labels to internal scatterer parameter names.
        """
        _axis_mapping = {}
        for tab in self.other_tabs:
            _axis_mapping.update(tab.component.mapping)

        return _axis_mapping
