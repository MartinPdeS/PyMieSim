#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
from tkinter import ttk
import tkinter
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widgets import InputWidget
from PyMieSim.gui.widget_collection import WidgetCollection


class ScattererTab(BaseTab):
    """
    A GUI tab for configuring scatterer parameters within PyMieSim experiments. This tab allows users
    to choose between different scatterer types (Sphere, Cylinder, CoreShell) and set relevant
    parameters like dimensions and refractive indices.

    Attributes:
        master (ttk.Notebook): The notebook widget this tab is part of.
        label (str): The label for the tab.
        source_tab (BaseTab): Reference to the source tab for source component configurations.
    """

    def __init__(self, master: ttk.Notebook, label: str, source_tab: BaseTab) -> NoReturn:
        """
        Initializes the ScattererTab with user interface elements for scatterer configuration.

        Args:
            master (ttk.Notebook): The notebook widget this tab will be part of.
            label (str): The label for the tab.
            source_tab (BaseTab): Reference to the source tab for accessing source component settings.
        """
        self.type_button = tkinter.StringVar(value='Sphere')
        self.source_tab = source_tab

        super().__init__(master, label=label)
        self._setup_combobox()
        self.setup_widgets()

    def setup_widgets(self) -> NoReturn:
        """
        Configures the GUI elements for the Scatterer tab based on the selected scatterer type.
        """
        match self.type_button.get():
            case 'Sphere':
                self.setup_sphere_widgets()
            case 'Cylinder':
                self.setup_cylinder_widgets()
            case 'CoreShell':
                self.setup_coreshell_widgets()
            case _:
                raise ValueError('Scatterer type not valid')

    def _setup_combobox(self) -> NoReturn:
        """
        Sets up a combobox for selecting the type of scatterer. It provides options for Sphere, Cylinder,
        or CoreShell configurations.
        """
        combobox = tkinter.ttk.Combobox(
            self.frame,
            textvariable=self.type_button,
            values=['Sphere', 'Cylinder', 'CoreShell'],
            state="readonly"
        )

        combobox.grid(row=0, column=0)
        combobox.bind("<<ComboboxSelected>>", self.on_type_change)

    def on_type_change(self, event=None) -> NoReturn:
        """
        Handles the event triggered by a change in the selected scatterer type. It updates the UI to match
        the selected scatterer configuration.

        Args:
            event: The event that triggered this method (default is None).
        """
        scatterer_type = self.type_button.get().lower()
        setup_method = getattr(self, f"setup_{scatterer_type}_widgets", None)
        self.widget_collection.clear_widgets()
        if callable(setup_method):
            setup_method()
        else:
            raise ValueError(f"Unsupported scatterer type: {scatterer_type}")

    def setup_sphere_widgets(self) -> NoReturn:
        """
        Sets up the configuration widgets for a Sphere scatterer and initializes the component.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(
            InputWidget(default_value='500', label='Diameter [nm]', component_label='diameter', multiplicative_factor=1e-9, dtype=float),
            InputWidget(default_value='1.4', label='Refractive Index', component_label='index', dtype=float),
            InputWidget(default_value='1.0', label='Medium Refractive Index', component_label='medium_index', dtype=float)
        )
        self.widget_collection.setup_widgets()
        self.setup_sphere_component()

    def setup_cylinder_widgets(self) -> NoReturn:
        """
        Sets up the configuration widgets for a Cylinder scatterer and initializes the component.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(
            InputWidget(default_value='1000', label='Diameter [nm]', component_label='diameter', multiplicative_factor=1e-9, dtype=float),
            InputWidget(default_value='1.4', label='Refractive Index', component_label='index', dtype=complex),
            InputWidget(default_value='1.0', label='Medium Refractive Index', component_label='medium_index', dtype=float)
        )
        self.widget_collection.setup_widgets()
        self.setup_cylinder_component()

    def setup_coreshell_widgets(self) -> NoReturn:
        """
        Sets up the configuration widgets for a CoreShell scatterer and initializes the component.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(
            InputWidget(default_value='1000', label='Core Diameter [nm]', component_label='core_diameter', multiplicative_factor=1e-9, dtype=float),
            InputWidget(default_value='200', label='Shell Width [nm]', component_label='shell_width', multiplicative_factor=1e-9, dtype=float),
            InputWidget(default_value='1.4', label='Core Refractive Index', component_label='core_index', dtype=complex),
            InputWidget(default_value='1.4', label='Shell Refractive Index', component_label='shell_index', dtype=complex),
            InputWidget(default_value='1.0', label='Medium Refractive Index', component_label='medium_index', dtype=float)
        )
        self.widget_collection.setup_widgets()
        self.setup_coreshell_component()

    def setup_component(self, event=None) -> NoReturn:
        """
        Handles the event triggered by a change in the selected scatterer type. It updates the UI to match
        the selected scatterer configuration.

        Args:
            event: The event that triggered this method (default is None).
        """
        scatterer_type = self.type_button.get().lower()
        self.widget_collection.update()
        setup_method = getattr(self, f"setup_{scatterer_type}_component", None)
        if callable(setup_method):
            setup_method()
        else:
            raise ValueError(f"Unsupported scatterer type: {scatterer_type}")

    def setup_sphere_component(self) -> NoReturn:
        self.component = scatterer.Sphere(**self.widget_collection.to_component_dict(), source=self.source_tab.component)

    def setup_cylinder_component(self) -> NoReturn:
        self.component = scatterer.Cylinder(**self.widget_collection.to_component_dict(), source=self.source_tab.component)

    def setup_coreshell_component(self) -> NoReturn:
        self.component = scatterer.CoreShell(**self.widget_collection.to_component_dict(), source=self.source_tab.component)

# -
