#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
from PyMieSim.experiment.source import Gaussian
from PyMieSim.gui.base_tab import BaseTab, Widget, WidgetCollection


class SourceTab(BaseTab):
    """
    The Source tab for configuring the light source in the simulation.

    Inherits from Tab and adds specific controls for configuring the source, such as
    wavelength, polarization, optical power, and numerical aperture.
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes the Source tab with predefined variables for the simulation's source configuration.
        """

        wavelength_widget = Widget(
            default_value='400:950:200',
            label='Wavelength [nm]',
            component_label='wavelength',
            multiplicative_factor=1e-9
        )

        polarization_widget = Widget(
            default_value='0',
            label='Polarization angle [degree]',
            component_label='polarization_value',
        )

        power_widget = Widget(
            default_value='1.0',
            label='Optical Power [mW]',
            component_label='optical_power',
            multiplicative_factor=1e-3
        )

        na_widget = Widget(
            default_value='0.2',
            label='Numerical Aperture (NA)',
            component_label='NA',
        )

        self.variables = WidgetCollection(wavelength_widget, polarization_widget, power_widget, na_widget)

        super().__init__(*args, **kwargs)

    def setup_tab(self) -> NoReturn:
        """
        Sets up the GUI elements for the Source tab, including labels and entry fields for each source parameter.
        """
        self.update_user_input()

        self.variables.setup_widgets(frame=self.frame)

        self.setup_component()

    def setup_component(self) -> NoReturn:
        self.update_user_input()

        self.component = Gaussian(**self.variables.to_component_dict())

        self.mapping = dict(
            wavelength=self.component.wavelength,
            polarization=self.component.polarization_value
        )
