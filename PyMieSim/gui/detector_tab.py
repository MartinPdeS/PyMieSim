#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
from tkinter import ttk, StringVar
import tkinter
from PyMieSim.experiment.detector import Photodiode, CoherentMode
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widgets import InputWidget, RadioButtonWidget
from PyMieSim.gui.widget_collection import WidgetCollection
from pydantic.dataclasses import dataclass
from pydantic import ConfigDict


@dataclass(kw_only=True, config=ConfigDict(arbitrary_types_allowed=True))
class DetectorTab(BaseTab):
    """
    A GUI tab for configuring the detector parameters for simulations in PyMieSim.

    Allows for the setup of detector characteristics such as numerical aperture (NA),
    angular offsets (gamma and phi), and polarization filter angle, facilitating
    the detailed configuration of the simulation's detector component.

    Attributes:
        variables (WidgetCollection): A collection of widgets for detector configuration.
        x_axis (tkinter.StringVar): empty.
        STD_axis (tkinter.StringVar): empty.

    Inherited attributes:
        notebook (ttk.Notebook): The notebook widget this tab is part of.
        label (str): The label for the tab.
        frame (ttk.Frame): The frame serving as the container for the tab's contents.
        main_window: Reference to the main window of the application, if applicable.
    """
    x_axis: tkinter.StringVar
    STD_axis: tkinter.StringVar

    def __post_init__(self) -> None:
        """
        Calls for BaseTab's post initialisation, and initializes the SourceTab with UI components for source configuration
        """
        super().__post_init__()
        self.type_button = StringVar(value='Photodiode')
        self.setup_type_combobox()
        self.setup_widgets()

    def setup_type_combobox(self) -> None:
        """
        Create and configure a combobox to select the type of detector, binding it to update UI on change.
        """
        self.type_widget = ttk.Combobox(
            self.frame,
            textvariable=self.type_button,
            values=['Photodiode', 'CoherentMode'],
            state="readonly"
        )
        self.type_widget.grid(row=0, column=0)
        self.type_widget.bind("<<ComboboxSelected>>", self.on_type_change)

    def on_type_change(self, event=None) -> NoReturn:
        """
        Handles the event triggered by a change in the selected scatterer type. It updates the UI to match
        the selected scatterer configuration.

        Args:
            event: The event that triggered this method (default is None).
        """
        detector_type = self.type_widget.get().lower()
        setup_method = getattr(self, f"setup_{detector_type}_widgets", None)
        self.widget_collection.clear_widgets()
        if callable(setup_method):
            setup_method()
        else:
            raise ValueError(f"Unsupported detector type: {detector_type}")

    def setup_widgets(self) -> NoReturn:
        """
        Configures the GUI elements for the Scatterer tab based on the selected scatterer type.
        """
        detector_type = self.type_widget.get()

        match detector_type:
            case 'Photodiode':
                self.setup_photodiode_widgets()
            case 'CoherentMode':
                self.setup_coherentmode_widgets()
            case _:
                raise ValueError(f'Detector type not valid: {detector_type}')

    def setup_photodiode_widgets(self) -> None:
        """
        Setup widgets specific to configuring a Photodiode detector.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(
            InputWidget(default_value='0.2, 0.3, 0.4', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Numerical aperture (NA)', component_label='NA', dtype=float),
            InputWidget(default_value='0', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Gamma [degree]', component_label='gamma_offset', dtype=float),
            InputWidget(default_value='0:360:200', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Phi [degree]', component_label='phi_offset', dtype=float),
            InputWidget(default_value='None', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Polarization filter [degree]', component_label='polarization_filter', dtype=float),
            InputWidget(default_value='500', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Sampling', component_label='sampling', dtype=int)
        )

        self.widget_collection.setup_widgets(row_start=1)
        self.setup_photodiode_component()

    def setup_coherentmode_widgets(self) -> None:
        """
        Setup widgets specific to configuring a Coherent Mode detector.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(
            RadioButtonWidget(option_text=['Point', 'Mean'], options_values=[False, True], component_label='mean_coupling', label='Mean coupling'),
            InputWidget(default_value='0', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Polarization filter [degree]', component_label='polarization_filter', dtype=float),
            InputWidget(default_value='0', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Gamma [degree]', component_label='gamma_offset', dtype=float),
            InputWidget(default_value='180:-180:200', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Phi [degree]', component_label='phi_offset', dtype=float),
            InputWidget(default_value='0.2, 0.3, 0.4', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Numerical aperture (NA)', component_label='NA', dtype=float),
            InputWidget(default_value='LP01', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Mode field', component_label='mode_number', dtype=str),
            InputWidget(default_value='0', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Field rotation [degree]', component_label='rotation', dtype=float),
            InputWidget(default_value='500', x_axis=self.x_axis, STD_axis=self.STD_axis, label='Sampling', component_label='sampling', dtype=int)
        )

        self.widget_collection.setup_widgets(row_start=1)
        self.setup_coherentmode_component()

    def setup_component(self, event=None) -> NoReturn:
        """
        Handles the event triggered by a change in the selected scatterer type. It updates the UI to match
        the selected scatterer configuration.

        Args:
            event: The event that triggered this method (default is None).
        """
        detector_type = self.type_button.get().lower()
        self.widget_collection.update()
        setup_method = getattr(self, f"setup_{detector_type}_component", None)
        if callable(setup_method):
            setup_method()
        else:
            raise ValueError(f"Unsupported scatterer type: {detector_type}")

    def setup_photodiode_component(self) -> NoReturn:
        kwargs = self.widget_collection.to_component_dict()

        self.component_dict = kwargs

        self.component = Photodiode(**kwargs)

    def setup_coherentmode_component(self) -> NoReturn:
        kwargs = self.widget_collection.to_component_dict()

        self.component_dict = kwargs

        self.component = CoherentMode(**kwargs)

# -
