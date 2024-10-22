#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
from pytest import raises
import tkinter
from PyMieSim.gui.main_window import PyMieSimGUI
from PyMieSim.gui.singleton import datashelf
from unittest.mock import patch


def set_up_gui(foo):
    """
    This is a decorator that will set up the gui, run the function and destroy the gui
    """
    def set_up():
        root = tkinter.Tk()
        root.geometry("750x600")
        PyMieSimGUI(root)

        foo()

        root.destroy()

    return set_up


@patch('tkinter.messagebox.showerror')
def radio_button_invoke(mock, widgets: list) -> None:
    for widget in widgets:
        if widget.can_be_axis:
            # Defining the radiobuttons
            radio_button_x_axis = widget.tk_radio_button_1
            radio_button_STD_axis = widget.tk_radio_button_2

            # The first part of the loop makes sure the buttons work individually
            # Checks if the x-axis radiobuttons work
            radio_button_x_axis.invoke()
            assert radio_button_x_axis['value'] == datashelf.x_axis_label_widget.get(), f"x-axis selection for the {radio_button_x_axis['value']} radio button did not work"
            datashelf.x_axis_label_widget.set(None)

            # Checks if the std-axis radiobuttons work
            radio_button_STD_axis.invoke()
            assert radio_button_STD_axis['value'] == datashelf.STD_axis_label_widget.get(), f"std-axis selection for the {radio_button_STD_axis['value']} radio button did not work"
            datashelf.x_axis_label_widget.set(None)

            # The second part of the loop checks if the correct ValueError gets raised if both selected axis are the same
            radio_button_x_axis.invoke()
            radio_button_STD_axis.invoke()
            with raises(ValueError):
                datashelf.control_tab.calculate_plot()


"""
The following three tests are checking if pressing the
radio buttons of the GUI is possible and if the variables self.STD_axis_label_widget and self.STD_axis_label_widget of the PyMieSimGUI class are updated upon clicking.
"""


@set_up_gui
def test_source_widgets() -> None:
    widgets = datashelf.source_tab.widget_collection.widgets
    radio_button_invoke(widgets=widgets)


@set_up_gui
def test_scatterer_widgets() -> None:
    for tab in datashelf.scatterer_tab.combobox_widget_collection.combobox_widget.tk_widget['values']:
        datashelf.scatterer_tab.combobox_widget_collection.combobox_widget.tk_widget.set(tab)
        datashelf.scatterer_tab.on_type_change()
        widgets = datashelf.scatterer_tab.widget_collection.widgets
        radio_button_invoke(widgets=widgets)


@set_up_gui
def test_detector_widgets() -> None:
    for tab in datashelf.detector_tab.combobox_widget_collection.combobox_widget.tk_widget['values']:
        datashelf.detector_tab.combobox_widget_collection.combobox_widget.tk_widget.set(tab)
        datashelf.detector_tab.on_type_change()
        widgets = datashelf.detector_tab.widget_collection.widgets
        radio_button_invoke(widgets=widgets)


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
