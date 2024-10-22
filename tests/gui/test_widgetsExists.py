#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import PyMieSim
from PyMieSim.gui.main_window import PyMieSimGUI
from PyMieSim.gui.singleton import datashelf
import tkinter


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


@set_up_gui
def test_input_widgets_exist():
    """
    This function checks if the number of widgets in each tab matches the expected number of widgets.
    """

    assert len(datashelf.source_tab.widget_collection.widgets) == 4, 'Missing widgets in the source_tab'

    for tab, widget_count in zip(datashelf.scatterer_tab.combobox_widget_collection.combobox_widget.tk_widget['values'], [3, 3, 5]):
        datashelf.scatterer_tab.combobox_widget_collection.combobox_widget.tk_widget.set(tab)
        datashelf.scatterer_tab.on_type_change()
        assert len(datashelf.scatterer_tab.widget_collection.widgets) == widget_count, f'Missing widget in the source/{tab} tab'

    for tab, widget_count in zip(datashelf.detector_tab.combobox_widget_collection.combobox_widget.tk_widget['values'], [5, 8]):
        datashelf.detector_tab.combobox_widget_collection.combobox_widget.tk_widget.set(tab)
        datashelf.detector_tab.on_type_change()
        assert len(datashelf.detector_tab.widget_collection.widgets) == widget_count, f'Missing widget in the detector/{tab} tab'


@set_up_gui
def test_control_button_exist():
    """
    This test ensures that an instance of all the following classes has been created by PyMieSimGUI.
    """
    assert datashelf.control_tab.calculate_button.__class__ == PyMieSim.gui.widgets.ControlWidget
    assert datashelf.control_tab.save_button.__class__ == PyMieSim.gui.widgets.ControlWidget
    assert datashelf.control_tab.export_button.__class__ == PyMieSim.gui.widgets.ControlWidget
    assert datashelf.control_tab.reset_std_button.__class__ == PyMieSim.gui.widgets.ControlWidget


if __name__ == "__main__":
    pytest.main(["-W error", __file__])
