#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn

from PyMieSim.experiment.source import Gaussian
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widget_collection import WidgetCollection
from PyMieSim.gui.singleton import datashelf

from pydantic.dataclasses import dataclass
from pydantic import ConfigDict
from tkinter import ttk


@dataclass(config=ConfigDict(arbitrary_types_allowed=True, kw_only=True))
class SourceTab(BaseTab):
    """
    A GUI tab for configuring the light source parameters for simulations in PyMieSim.

    This class provides a user interface for setting up the light source by specifying
    parameters such as wavelength, polarization, optical power, and numerical aperture (NA).
    User inputs are used to configure a Gaussian light source in the simulation.

    Attributes:
        variables (WidgetCollection): A collection of widgets for source configuration.

    Inherited attributes:
        notebook (ttk.Notebook): The notebook widget this tab is part of.
        label (str): The label for the tab.
        main_window: Reference to the main window of the application, if applicable.
    """
    notebook: ttk.Notebook
    label: str
    main_window = None

    def __post_init__(self):
        """
        Calls for BaseTab's post initialisation, and initializes the SourceTab with UI components for source configuration
        """
        super().__init__()
        self.setup_widgets()

    def setup_widgets(self) -> NoReturn:
        """
        Configures the GUI elements for the Source tab.

        This method sets up labels and entry fields for each source parameter, facilitating
        user interaction for the configuration of the light source.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(tab='source_tab', component='Gaussian')

        self.widget_collection.setup_widgets()
        self.setup_component()

    def setup_component(self) -> NoReturn:
        """
        Initializes the Gaussian source component based on user input.

        This method reads input values from the UI widgets and uses them to configure
        a Gaussian source component for the simulation.
        """
        self.widget_collection.update()
        kwargs = self.widget_collection.to_component_dict()

        kwargs["optical_power"] = kwargs["optical_power"][0]

        self.component = Gaussian(**kwargs)
        datashelf.source_component = self.component


# -
