#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn, Dict
import tkinter

from PyMieSim.gui.widgets import BaseWidget
from PyMieSim.gui.widget_dictonary import widget_dock


class WidgetCollection:
    """
    A collection class for managing multiple Widget instances.

    This class facilitates the grouping of Widget instances, allowing for collective operations like
    updating widget values, clearing non-permanent widgets, and setting up widgets within a tkinter frame.

    Attributes:
        widgets (tuple[Widget, ...]): A tuple of Widget instances included in the collection.
    """

    def __init__(self, frame, *widgets: BaseWidget) -> None:
        """
        Initializes a new instance of WidgetCollection with a given set of Widget instances.

        Args:
            *widgets (Widget): Variable length Widget instances to be included in the collection.
        """
        self.frame = frame
        self.row_start = 0

    def add_widgets(self, tab: str, component: str):
        self.widgets = widget_dock[tab][component]

        for widget in self.widgets:
            widget.frame = self.frame
            widget.initialize()

    def setup_combobox_widget(self, tab: str, component: str):
        self.combobox_widget = widget_dock[tab][component][0]
        self.combobox_widget.frame = self.frame
        self.combobox_widget.setup()

    def setup_control_widget(self, config_dict, tab: str = 'control_tab'):
        column_counter = 0

        self.widgets = widget_dock[tab]

        for self.widget, self.component_label in zip(self.widgets, config_dict.keys()):
            self.widget.component_label = self.component_label
            self.widget.frame = self.frame
            self.widget.command = config_dict[self.component_label]
            self.widget.setup(column=column_counter)
            column_counter += 1

    def to_component_dict(self) -> Dict[str, float | str]:
        """
        Creates a dictionary mapping component labels to their respective widget values.

        Returns:
            dict[str, float | str]: A dictionary where keys are component labels and values are widget values.
        """
        return {widget.component_label: widget.get_value() for widget in self.widgets}

    def __getitem__(self, component_label: str) -> BaseWidget:
        """
        Allows direct access to a Widget instance in the collection by its component label.

        Args:
            component_label (str): The component label of the desired Widget.

        Returns:
            Widget: The Widget instance with the matching component label, if found.
        """
        return next((widget for widget in self.widgets if widget.component_label == component_label), None)

    def clear_widgets(self) -> NoReturn:
        """
        Clears all non-permanent widgets from the tkinter frame.
        """
        for widget in self.widgets:
            if not widget.is_permanent:
                widget.destroy()
                widget.destroy()

    def update(self) -> NoReturn:
        """
        Updates the value of all widgets in the collection.
        """
        for widget in self.widgets:
            widget.update()

    def title_bar(self, title_bar: bool = False):
        """
        Sets up a title bar on each page of the source, detector and scatterer tabs

        Args:
            title_bar: A boolean indiquating if the tab requires a. titble_bar
        """
        if title_bar:
            self.tk_label = tkinter.Label(self.frame, text="Variable")
            self.tk_label.grid(row=1, column=0, sticky="W", pady=2)
            self.tk_widget_title = tkinter.Label(self.frame, text="Values")
            self.tk_widget_title.grid(row=1, column=1, sticky="W", pady=2)
            self.x_axis_title = tkinter.Label(self.frame, text="x-axis")
            self.x_axis_title.grid(row=1, column=2, sticky="W", pady=2)
            self.STD_axis_title = tkinter.Label(self.frame, text="STD-axis")
            self.STD_axis_title.grid(row=1, column=3, sticky="W", pady=2)

    def setup_widgets(self, row_start: int = 1, title_bar=True) -> NoReturn:
        """
        Sets up and packs the widgets within a specified tkinter frame.

        Args:
            frame (tkinter.Frame): The tkinter frame where widgets are to be packed.
        """
        self.title_bar(title_bar)

        for row, widget in enumerate(self.widgets):
            widget.setup(row=row + row_start)

    def __repr__(self) -> str:
        """
        Returns a string representation of the WidgetCollection instance, listing all widget labels.

        Returns:
            str: A space-separated string of widget labels.
        """
        return " ".join(str(widget) for widget in self.widgets)

# -
