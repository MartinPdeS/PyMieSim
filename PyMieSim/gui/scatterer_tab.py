#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn, Dict, Union
import tkinter
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.utils import Widget, WidgetCollection


class ScattererTab(BaseTab):
    """
    GUI tab for configuring the scatterer parameters in a PyMieSim experiment.

    Allows users to select the scatterer type (Sphere, Cylinder, CoreShell) and specify relevant parameters,
    including dimensions and refractive indices. This class manages user input and updates the scatterer configuration
    accordingly.

    Attributes:
        type_button (tkinter.StringVar): Holds the selected scatterer type.
        self.sphere_widget_collection (WidgetCollection): Widgets for configuring a Sphere scatterer.
        self.cylinder_widget_collection (WidgetCollection): Widgets for configuring a Cylinder scatterer.
        self.coreshell_widget_collection (WidgetCollection): Widgets for configuring a CoreShell scatterer.
        source_tab: Reference to the source tab for accessing the source component.
    """

    def __init__(self, *args, source_tab, **kwargs):
        """
        Initializes the Scatterer tab with default scatterer configuration widgets and binds events.

        Parameters:
            source_tab: The tab that contains the source component settings.
            *args, **kwargs: Additional arguments and keyword arguments for the BaseTab.
        """
        self.type_button = tkinter.StringVar(value='Sphere')

        self.sphere_widget_collection = WidgetCollection(
            Widget(default_value='1000', label='Diameter [nm]', component_label='diameter', multiplicative_factor=1e-9),
            Widget(default_value='1.4', label='Refractive Index', component_label='index'),
            Widget(default_value='1.0', label='Medium Refractive Index', component_label='n_medium')
        )

        self.cylinder_widget_collection = WidgetCollection(
            Widget(default_value='1000', label='Diameter [nm]', component_label='diameter', multiplicative_factor=1e-9),
            Widget(default_value='1.4', label='Refractive Index', component_label='index'),
            Widget(default_value='1.0', label='Medium Refractive Index', component_label='n_medium')
        )

        self.coreshell_widget_collection = WidgetCollection(
            Widget(default_value='1000', label='Core Diameter [nm]', component_label='core_diameter', multiplicative_factor=1e-9),
            Widget(default_value='1000', label='Shell Width [nm]', component_label='shell_width', multiplicative_factor=1e-9),
            Widget(default_value='1.4', label='Core Refractive Index', component_label='core_index'),
            Widget(default_value='1.4', label='Shell Refractive Index', component_label='shell_index'),
            Widget(default_value='1.0', label='Medium Refractive Index', component_label='n_medium')
        )

        self.source_tab = source_tab

        super().__init__(*args, **kwargs)

        combobox = tkinter.ttk.Combobox(
            self.frame,
            textvariable=self.type_button,
            values=['Sphere', 'Cylinder', 'CoreShell'],
            state="readonly"
        )

        combobox.pack(side=tkinter.TOP)

        combobox.bind("<<ComboboxSelected>>", self.on_type_change)

    def on_type_change(self, event=None) -> NoReturn:
        """
        Handles the scatterer type change event, updates the UI widgets and scatterer configuration accordingly.
        """
        self.widget_collection.clear_widgets()

        self.setup_tab()

        self.main_window.axis_tab.clear_button()

        self.main_window.axis_tab.setup_tab()

    def get_selected_type(self) -> str:
        """
        Retrieves the currently selected scatterer type.

        Returns:
            The selected scatterer type as a string.
        """
        return self.type_button.get()

    def get_variables(self) -> WidgetCollection:
        """
        Retrieves the widget collection corresponding to the selected scatterer type.

        Returns:
            A WidgetCollection instance for the selected scatterer type.
        """
        match self.get_selected_type():
            case 'Sphere':
                return self.sphere_widget_collection
            case 'Cylinder':
                return self.cylinder_widget_collection
            case 'CoreShell':
                return self.coreshell_widget_collection
            case _:
                raise ValueError('Scatterer type not valid')

    def setup_tab(self) -> NoReturn:
        """
        Configures the GUI elements for the Scatterer tab based on the selected scatterer type.
        """
        self.widget_collection = self.get_variables()

        self.widget_collection.setup_widgets(frame=self.frame)

        self.setup_component()

    def setup_sphere_component(self) -> NoReturn:
        """Configures a sphere scatterer component with the selected parameters."""
        self.component = scatterer.Sphere(
            **self.widget_collection.to_component_dict(),
            source_set=self.source_tab.component
        )

        self.mapping = dict(
            diameter=self.component.diameter,
            index=self.component.index,
            n_medium=self.component.n_medium,
        )

    def setup_cylinder_component(self) -> NoReturn:
        """Configures a cylinder scatterer component with the selected parameters."""
        self.component = scatterer.Cylinder(
            **self.widget_collection.to_component_dict(),
            source_set=self.source_tab.component
        )

        self.mapping = dict(
            diameter=self.component.diameter,
            index=self.component.index,
            n_medium=self.component.n_medium,
        )

    def setup_coreshell_component(self) -> NoReturn:
        """Configures a core-shell scatterer component with the selected parameters."""
        self.component = scatterer.CoreShell(
            **self.widget_collection.to_component_dict(),
            source_set=self.source_tab.component
        )

        self.mapping = dict(
            core_diameter=self.component.core_diameter,
            shell_width=self.component.shell_width,
            core_index=self.component.core_index,
            shell_index=self.component.shell_index,
            n_medium=self.component.n_medium,
        )

    def setup_component(self) -> NoReturn:
        self.update_user_input()

        match self.get_selected_type():
            case 'Sphere':
                self.setup_sphere_component()
            case 'Cylinder':
                self.setup_cylinder_component()
            case 'CoreShell':
                self.setup_coreshell_component()
            case _:
                raise ValueError('Scatterer type not valid')

# -
