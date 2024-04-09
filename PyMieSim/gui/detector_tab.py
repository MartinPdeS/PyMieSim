#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import tkinter
from PyMieSim.experiment.detector import Photodiode, LPMode
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.utils import InputWidget, WidgetCollection


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
        self.type_button = tkinter.StringVar(value='Photodiode')

        super().__init__(*args, **kwargs)

        self._setup_combobox()

    def _setup_combobox(self):
        """
        Sets up the combobox for selecting the scatterer type.
        """
        combobox = tkinter.ttk.Combobox(
            self.frame,
            textvariable=self.type_button,
            values=['Photodiode', 'LPMode'],
            state="readonly"
        )

        combobox.pack(side=tkinter.TOP)
        combobox.bind("<<ComboboxSelected>>", self.on_type_change)

    def on_type_change(self, event=None) -> NoReturn:
        """
        Handles the scatterer type change event, updates the UI widgets and scatterer configuration accordingly.
        """
        self.widget_collection.clear_widgets()

        self.setup()

    def get_selected_type(self) -> str:
        """
        Retrieves the currently selected scatterer type.

        Returns:
            The selected scatterer type as a string.
        """
        return self.type_button.get()

    def setup(self) -> NoReturn:
        """
        Configures the GUI elements for the Detector tab.

        Sets up labels, entry fields, and other UI components based on the detector
        configuration options, enabling user interaction for configuring the detector.
        """
        match self.get_selected_type().lower():
            case 'photodiode':
                self._setup_photodiode_widgets()
                self._setup_photodiode_component()
            case 'lpmode':
                self._setup_lpmode_widgets()
                self._setup_lpmode_component()
            case _:
                raise ValueError('Detector type not valid')

    def _setup_photodiode_widgets(self):
        self.widget_collection = WidgetCollection(
            InputWidget(default_value='0.2, 0.3, 0.4', label='Numerical aperture (NA)', component_label='NA'),
            InputWidget(default_value='0', label='Gamma [degree]', component_label='gamma_offset'),
            InputWidget(default_value='0:360:200', label='Phi [degree]', component_label='phi_offset'),
            InputWidget(default_value='None', label='Polarization filter [degree]', component_label='polarization_filter')
        )

        self.widget_collection.setup_widgets(frame=self.frame)

    def _setup_lpmode_widgets(self):
        self.widget_collection = WidgetCollection(
            # InputWidget(default_value='point', label='Coupling mode', component_label='coupling_mode', to_float=False),
            InputWidget(default_value='0.', label='Polarization filter [degree]', component_label='polarization_filter'),
            InputWidget(default_value='0', label='Gamma [degree]', component_label='gamma_offset'),
            InputWidget(default_value='180:-180:200', label='Phi [degree]', component_label='phi_offset'),
            InputWidget(default_value='0.2, 0.3, 0.4', label='Numerical aperture (NA)', component_label='NA'),
            InputWidget(default_value='LP01', label='Mode field', component_label='mode_number', to_float=False),
        )

        self.widget_collection.setup_widgets(frame=self.frame)

    def _setup_photodiode_component(self):
        self.component = Photodiode(**self.widget_collection.to_component_dict(), sampling=500)

        self.mapping = {
            'NA': self.component.NA,
            'gamma': self.component.gamma_offset,
            'phi': self.component.phi_offset,
            'polarization_filter': self.component.polarization_filter
        }

    def _setup_lpmode_component(self):
        self.component = LPMode(**self.widget_collection.to_component_dict(), sampling=500)

        self.mapping = {
            'NA': self.component.NA,
            'gamma': self.component.gamma_offset,
            'phi': self.component.phi_offset,
            'polarization_filter': self.component.polarization_filter
        }

    def setup_component(self) -> NoReturn:
        """
        Initializes the Photodiode detector component based on user input.

        This method reads input values from the UI widgets and uses them to configure
        a Photodiode detector component for the simulation, including setting the sampling
        parameter explicitly.
        """
        self.update_user_input()

        selected_type = self.get_selected_type()

        match selected_type.lower():
            case 'photodiode':
                self._setup_photodiode_component()
            case 'lpmode':
                self._setup_lpmode_component()

# -
