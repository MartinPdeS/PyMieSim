#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn

from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widget_collection import WidgetCollection

from pydantic.dataclasses import dataclass
from pydantic import ConfigDict
from tkinter import ttk


@dataclass(config=ConfigDict(arbitrary_types_allowed=True, kw_only=True))
class AxisTab(BaseTab):
    """
    A GUI tab for selecting the y axis of the PyMieSim experiment, as well as configuring its x and std axis.

    Attributes:
        other_tabs (list[BaseTab]): The source, scatterer and detector tabs, in which the x and std axis are selected

    Inherited attributes:
        notebook (ttk.Notebook): The notebook widget this tab is part of.
        label (str): The label for the tab.
        frame (ttk.Frame): The frame serving as the container for the tab's contents.
        main_window: Reference to the main window of the application, if applicable.
    """
    notebook: ttk.Notebook
    label: str
    main_window = None

    def __post_init__(self) -> NoReturn:
        """
         Calls for BaseTab's post initialisation, and initializes the Axis Configuration tab with references to other tabs to gather possible axis choices.
        """
        super().__init__()
        self.setup()

    def setup(self) -> NoReturn:
        """
        Sets up the UI elements for axis configuration, including comboboxes for selecting
        variables for the x-axis, y-axis, and an optional standard deviation (STD) axis.
        """

        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(tab='axis_tab', component='y_axis')

        self.widget_collection.setup_widgets(title_bar=False)

    def update(self) -> NoReturn:
        self.widget_collection.clear_widgets()

        self.setup()

# -
