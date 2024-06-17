#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn, Dict
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widgets import ComBoxWidget
from PyMieSim.gui.widget_collection import WidgetCollection
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict


@dataclass(kw_only=True, config=ConfigDict(arbitrary_types_allowed=True))
class AxisTab(BaseTab):
    """
    A GUI tab for selecting the y axis of the PyMieSim experiment, as well as configuring its x and std axis.

    Attributes:
        other_tabs (list[BaseTab]): The source, scatterer and detector tabs, in which the x and std axis are selected
        measure_map: The list of the different measures (i.e. y axis) possible to choose from

    Inherited attributes:
        notebook (ttk.Notebook): The notebook widget this tab is part of.
        label (str): The label for the tab.
        frame (ttk.Frame): The frame serving as the container for the tab's contents.
        main_window: Reference to the main window of the application, if applicable.
    """
    other_tabs: list[BaseTab]
    measure_map = scatterer.Sphere.available_measure_list  # Note: this list should be updated with the specific scatterer tab selected, instead of always using the sphere

    def __post_init__(self) -> NoReturn:
        """
         Calls for BaseTab's post initialisation, and initializes the Axis Configuration tab with references to other tabs to gather possible axis choices.
        """
        super().__post_init__()
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
            ComBoxWidget(label='y-axis', component_label='y_axis', options=self.y_axis_options, default_options=len(self.y_axis_options) - 1),
        )

        self.widget_collection.setup_widgets(title_bar=False)

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

# -
