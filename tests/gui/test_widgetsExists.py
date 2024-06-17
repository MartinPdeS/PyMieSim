#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyMieSim.gui.main_window import PyMieSimGUI
import tkinter as tk
from tkinter.ttk import Button


def set_up_gui(foo):
    """
    This is a decorator that will set up the gui, run the function and destroy the gui
    """
    def set_up():
        root = tk.Tk()
        root.geometry("750x600")
        gui = PyMieSimGUI(root)

        foo(gui=gui)

        root.destroy()

    return set_up


@set_up_gui
def test_input_widgets_exist(**kwargs):
    """
    This function checks if the number of widgets in each tab matches the expected number of widgets.
    """
    gui = kwargs['gui']

    assert len(gui.source_tab.widget_collection.widgets) == 4, 'Missing widgets in the source_tab'

    for tab, widget_count in zip(gui.scatterer_tab.type_widget['values'], [3, 3, 5]):
        gui.scatterer_tab.type_widget.set(tab)
        gui.scatterer_tab.on_type_change()
        assert len(gui.scatterer_tab.widget_collection.widgets) == widget_count, f'Missing widget in the source/{tab} tab'

    for tab, widget_count in zip(gui.detector_tab.type_widget['values'], [5, 8]):
        gui.detector_tab.type_widget.set(tab)
        gui.detector_tab.on_type_change()
        assert len(gui.detector_tab.widget_collection.widgets) == widget_count, f'Missing widget in the detector/{tab} tab'


@set_up_gui
def test_control_button_exist(**kwargs):
    """
    This test ensures that an instance of all the following classes has been created by PyMieSimGUI.
    """
    gui = kwargs['gui']
    assert gui.calculate_button.__class__ == Button
    assert gui.save_button.__class__ == Button
    assert gui.export_button.__class__ == Button
    assert gui.reset_std_button.__class__ == Button

# -
