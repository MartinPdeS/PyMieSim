#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import tkinter
from PyMieSim.experiment.source import Gaussian
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widgets import InputWidget
from PyMieSim.gui.widget_collection import WidgetCollection
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict


@dataclass(kw_only=True, config=ConfigDict(arbitrary_types_allowed=True))
class SourceTab(BaseTab):
    """
    A GUI tab for configuring the light source parameters for simulations in PyMieSim.

    This class provides a user interface for setting up the light source by specifying
    parameters such as wavelength, polarization, optical power, and numerical aperture (NA).
    User inputs are used to configure a Gaussian light source in the simulation.

    Attributes:
        variables (WidgetCollection): A collection of widgets for source configuration.
        x_axis (tkinter.StringVar): empty.
        STD_axis (tkinter.StringVar): empty.

    """
    x_axis: tkinter.StringVar
    STD_axis: tkinter.StringVar

    def __post_init__(self):
        """
        Calls for BaseTab's post initialisation, and initializes the SourceTab with UI components for source configuration
        """
        super().__post_init__()
        self.setup_widgets()

    def setup_widgets(self) -> NoReturn:
        """
        Configures the GUI elements for the Source tab.

        This method sets up labels and entry fields for each source parameter, facilitating
        user interaction for the configuration of the light source.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(
            InputWidget(default_value='1310', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Wavelength [nm]', component_label='wavelength', multiplicative_factor=1e-9, dtype=float),
            InputWidget(default_value='0', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Polarization angle [degree]', component_label='polarization', dtype=float),
            InputWidget(default_value='1.0', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Optical Power [mW] [fix]', component_label='optical_power',
                        multiplicative_factor=1e-3, can_be_axis=False, dtype=float),  # If can_be_axis is false, then will not put the yo widget!
            InputWidget(default_value='0.2', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Numerical Aperture (NA) [fix]', component_label='NA', can_be_axis=False, dtype=float),
        )

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

# -
