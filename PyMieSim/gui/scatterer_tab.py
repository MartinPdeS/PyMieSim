#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import tkinter
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.utils import InputWidget, WidgetCollection


class ScattererTab(BaseTab):
    """
    GUI tab for configuring scatterer parameters in PyMieSim experiments. It allows for the selection of
    scatterer type (Sphere, Cylinder, CoreShell) and specification of relevant parameters, including dimensions
    and refractive indices. The class manages user input to update scatterer configurations dynamically.

    Attributes:
        type_button (tk.StringVar): Variable to hold the selected scatterer type.
        widget_collections (dict): Dictionary of WidgetCollections for different scatterer configurations.
        source_tab: Reference to the source tab for accessing the source component configurations.
    """

    def __init__(self, *args, source_tab, **kwargs):
        """
        Initializes the ScattererTab with UI elements for scatterer configuration.

        Parameters:
            notebook (ttk.Notebook): The notebook widget this tab will be part of.
            label (str): The label for the tab.
            source_tab: The tab containing source component settings.
            **kwargs: Additional keyword arguments.
        """
        self.type_button = tkinter.StringVar(value='Sphere')

        self.source_tab = source_tab

        super().__init__(*args, **kwargs)

        self._setup_combobox()

    def _setup_combobox(self):
        """
        Sets up the combobox for selecting the scatterer type.
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
        Handles the scatterer type change event, updates the UI widgets and scatterer configuration accordingly.
        """
        self.widget_collection.clear_widgets()

        self.setup()

        self.main_window.axis_tab.clear_button()

        self.main_window.axis_tab.setup_tab()

    def get_selected_type(self) -> str:
        """
        Retrieves the currently selected scatterer type.

        Returns:
            The selected scatterer type as a string.
        """
        return self.type_button.get()

    def _setup_sphere_widgets(self) -> NoReturn:
        """Configures a sphere scatterer component with the selected parameters."""
        self.widget_collection = WidgetCollection(
            InputWidget(default_value='500', label='Diameter [nm]', component_label='diameter', multiplicative_factor=1e-9),
            InputWidget(default_value='1.4', label='Refractive Index', component_label='index'),
            InputWidget(default_value='1.0', label='Medium Refractive Index', component_label='n_medium')
        )

        self.widget_collection.setup_widgets(frame=self.frame)

    def setup_sphere_component(self) -> NoReturn:
        self.component = scatterer.Sphere(
            **self.widget_collection.to_component_dict(),
            source_set=self.source_tab.component
        )

        self.mapping = dict(
            diameter=self.component.diameter,
            index=self.component.index,
            n_medium=self.component.n_medium,
        )

    def _setup_cylinder_widgets(self) -> NoReturn:
        """Configures a cylinder scatterer component with the selected parameters."""
        self.widget_collection = WidgetCollection(
            InputWidget(default_value='1000', label='Diameter [nm]', component_label='diameter', multiplicative_factor=1e-9),
            InputWidget(default_value='1.4', label='Refractive Index', component_label='index'),
            InputWidget(default_value='1.0', label='Medium Refractive Index', component_label='n_medium')
        )

        self.widget_collection.setup_widgets(frame=self.frame)

    def setup_cylinder_component(self) -> NoReturn:
        self.component = scatterer.Cylinder(
            **self.widget_collection.to_component_dict(),
            source_set=self.source_tab.component
        )

        self.mapping = dict(
            diameter=self.component.diameter,
            index=self.component.index,
            n_medium=self.component.n_medium,
        )

    def _setup_coreshell_widgets(self) -> NoReturn:
        """Configures a core-shell scatterer component with the selected parameters."""
        self.widget_collection = WidgetCollection(
            InputWidget(default_value='1000', label='Core Diameter [nm]', component_label='core_diameter', multiplicative_factor=1e-9),
            InputWidget(default_value='1000', label='Shell Width [nm]', component_label='shell_width', multiplicative_factor=1e-9),
            InputWidget(default_value='1.4', label='Core Refractive Index', component_label='core_index'),
            InputWidget(default_value='1.4', label='Shell Refractive Index', component_label='shell_index'),
            InputWidget(default_value='1.0', label='Medium Refractive Index', component_label='n_medium')
        )

        self.widget_collection.setup_widgets(frame=self.frame)

    def setup_coreshell_component(self) -> NoReturn:
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

    def setup(self) -> NoReturn:
        """
        Configures the GUI elements for the Scatterer tab based on the selected scatterer type.
        """
        match self.get_selected_type():
            case 'Sphere':
                self._setup_sphere_widgets()
                self.setup_sphere_component()
            case 'Cylinder':
                self._setup_cylinder_widgets()
                self.setup_cylinder_component()
            case 'CoreShell':
                self._setup_coreshell_widgets()
                self.setup_coreshell_component()
            case _:
                raise ValueError('Scatterer type not valid')

    def setup_component(self) -> NoReturn:
        """
        Configures the GUI elements for the Scatterer tab based on the selected scatterer type.
        """
        self.update_user_input()

        selected_type = self.get_selected_type()

        match selected_type.lower():
            case 'sphere':
                self.setup_sphere_component()
            case 'cylinder':
                self.setup_cylinder_component()
            case 'coreshell':
                self.setup_coreshell_component()
            case _:
                raise ValueError('Scatterer type not valid')

# -
