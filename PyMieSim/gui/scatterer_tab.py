#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import tkinter

from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab
from PyMieSim.gui.widget_collection import WidgetCollection
from PyMieSim.gui.axis_tab import AxisTab
from PyMieSim.gui.singleton import datashelf

from pydantic.dataclasses import dataclass
from pydantic import ConfigDict
from tkinter import ttk


@dataclass(config=ConfigDict(arbitrary_types_allowed=True, kw_only=True))
class ScattererTab(BaseTab):
    """
    A GUI tab for configuring scatterer parameters within PyMieSim experiments. This tab allows users
    to choose between different scatterer types (Sphere, Cylinder, CoreShell) and set relevant
    parameters like dimensions and refractive indices.

    Attribute:
        axis_tab (AxisTab): Used to update the y axis selection in the axis tab when the scatterer changes

    Inherited attributes:
        notebook (ttk.Notebook): The notebook widget this tab is part of.
        label (str): The label for the tab.
        frame (ttk.Frame): The frame serving as the container for the tab's contents.
        main_window: Reference to the main window of the application, if applicable.
    """
    notebook: ttk.Notebook
    label: str
    main_window = None
    axis_tab: AxisTab

    def __post_init__(self) -> NoReturn:
        """
        Calls for BaseTab's post initialisation, and initializes the SourceTab with UI components for source configuration
        """
        self.type_button = tkinter.StringVar(value='Sphere')
        super().__init__()
        self._setup_combobox()
        self.setup_widgets()

    def setup_widgets(self) -> NoReturn:
        """
        Configures the GUI elements for the Scatterer tab based on the selected scatterer type.
        """
        match datashelf.scatterer_selection.get():
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
        self.combobox_widget_collection = WidgetCollection(frame=self.frame)

        self.combobox_widget_collection.setup_combobox_widget(tab='scatterer_tab', component='Combox')

        self.combobox_widget_collection.combobox_widget.tk_widget.bind("<<ComboboxSelected>>", self.on_type_change)

    def on_type_change(self, event=None) -> NoReturn:
        """
        Handles the event triggered by a change in the selected scatterer type. It updates the UI to match
        the selected scatterer configuration.

        Args:
            event: The event that triggered this method (default is None).
        """
        # Changing axis tab
        scatterer_type = datashelf.scatterer_selection.get()
        setup_method = getattr(self, f"setup_{scatterer_type.lower()}_widgets", None)

        self.widget_collection.clear_widgets()
        if callable(setup_method):
            setup_method()
        else:
            raise ValueError(f"Unsupported scatterer type: {scatterer_type}")

        # Changing y_axis options in accordance
        scatterer_type_class = getattr(scatterer, scatterer_type)
        datashelf.measure_map = getattr(scatterer_type_class, 'available_measure_list')
        self.axis_tab.update()

    def setup_sphere_widgets(self) -> NoReturn:
        """
        Sets up the configuration widgets for a Sphere scatterer and initializes the component.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(tab='scatterer_tab', component='Sphere')

        self.widget_collection.setup_widgets()
        self.setup_sphere_component()

    def setup_cylinder_widgets(self) -> NoReturn:
        """
        Sets up the configuration widgets for a Cylinder scatterer and initializes the component.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(tab='scatterer_tab', component='Cylinder')

        self.widget_collection.setup_widgets()
        self.setup_cylinder_component()

    def setup_coreshell_widgets(self) -> NoReturn:
        """
        Sets up the configuration widgets for a CoreShell scatterer and initializes the component.
        """
        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.add_widgets(tab='scatterer_tab', component='CoreShell')

        self.widget_collection.setup_widgets()
        self.setup_coreshell_component()

    def setup_component(self, event=None) -> NoReturn:
        """
        Handles the event triggered by a change in the selected scatterer type. It updates the UI to match
        the selected scatterer configuration.

        Args:
            event: The event that triggered this method (default is None).
        """
        scatterer_type = datashelf.scatterer_selection.get()
        self.widget_collection.update()
        setup_method = getattr(self, f"setup_{scatterer_type.lower()}_component", None)
        if callable(setup_method):
            setup_method()
        else:
            raise ValueError(f"Unsupported scatterer type: {scatterer_type}")

    def setup_sphere_component(self) -> NoReturn:
        kwargs = self.widget_collection.to_component_dict()
        self.component = scatterer.Sphere(**kwargs, source=datashelf.source_component)

    def setup_cylinder_component(self) -> NoReturn:
        kwargs = self.widget_collection.to_component_dict()
        self.component = scatterer.Cylinder(**kwargs, source=datashelf.source_component)

    def setup_coreshell_component(self) -> NoReturn:
        kwargs = self.widget_collection.to_component_dict()
        self.component = scatterer.CoreShell(**kwargs, source=datashelf.source_component)


# -
