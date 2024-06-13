#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn, Dict

from PyMieSim.gui.widgets import BaseWidget


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

    def add_widgets(self, *widgets) -> NoReturn:
        for widget in widgets:
            widget.frame = self.frame

        self.widgets = widgets

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

    def setup_widgets(self, row_start: int = 0) -> NoReturn:
        """
        Sets up and packs the widgets within a specified tkinter frame.

        Args:
            frame (tkinter.Frame): The tkinter frame where widgets are to be packed.
        """
        for row, widget in enumerate(self.widgets):
            widget.setup(row=row + row_start)

    def __repr__(self) -> str:
        """
        Returns a string representation of the WidgetCollection instance, listing all widget labels.

        Returns:
            str: A space-separated string of widget labels.
        """
        return " ".join(str(widget) for widget in self.widgets)
