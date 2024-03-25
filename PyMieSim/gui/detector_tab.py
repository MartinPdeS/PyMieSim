#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.utils import Widget, WidgetCollection


class DetectorTab(BaseTab):
    """
    A GUI tab for configuring the detector parameters for simulations in PyMieSim.

    Allows for the setup of detector characteristics such as numerical aperture (NA),
    angular offsets (gamma and phi), and polarization filter angle, facilitating
    the detailed configuration of the simulation's detector component.

    Attributes:
        variables (WidgetCollection): A collection of widgets for detector configuration.
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes the DetectorTab with UI components for detector configuration.

        Parameters:
            *args: Variable length argument list for BaseTab.
            **kwargs: Arbitrary keyword arguments for BaseTab.
        """
        self._init_widgets()

        super().__init__(*args, **kwargs)

    def _init_widgets(self):
        """Initializes widgets for detector configuration with default values and labels."""
        self.widget_collection = WidgetCollection(
            Widget(default_value='0.2, 0.4', label='Numerical aperture (NA)', component_label='NA'),
            Widget(default_value='0', label='Gamma [degree]', component_label='gamma_offset'),
            Widget(default_value='0', label='Phi [degree]', component_label='phi_offset'),
            Widget(default_value='0.', label='Polarization filter [degree]', component_label='polarization_filter')
        )

    def setup_tab(self) -> NoReturn:
        """
        Configures the GUI elements for the Detector tab.

        Sets up labels, entry fields, and other UI components based on the detector
        configuration options, enabling user interaction for configuring the detector.
        """
        self.widget_collection.setup_widgets(frame=self.frame)
        self.setup_component()

    def setup_component(self) -> NoReturn:
        """
        Initializes the Photodiode detector component based on user input.

        This method reads input values from the UI widgets and uses them to configure
        a Photodiode detector component for the simulation, including setting the sampling
        parameter explicitly.
        """
        self.update_user_input()
        self.component = Photodiode(**self.widget_collection.to_component_dict(), sampling=300)
        self._update_mapping()

    def _update_mapping(self):
        """
        Updates the mapping attribute with current component values.

        This private method extracts relevant parameters from the configured Photodiode
        detector component and updates the mapping attribute to reflect these values.
        """
        self.mapping = {
            'NA': self.component.NA,
            'gamma': self.component.gamma_offset,
            'phi': self.component.phi_offset,
            'polarization_filter': self.component.polarization_filter
        }

# -
