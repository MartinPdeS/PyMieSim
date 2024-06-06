from PyMieSim.gui.main_window import PyMieSimGUI
import tkinter as tk
from pytest import raises
from unittest.mock import patch


root = tk.Tk()
root.geometry("750x600")
gui = PyMieSimGUI(root)

"""
The following 3 tests are checking if pressing the radio buttons of the gui is 
possible and if the self.STD_axis_label_widget and self.STD_axis_label_widget variables of the
PyMieSimGUI class are updated on click
"""

def radio_button_invoke(widgets : list) -> None:
    """
    This function will take a list of widgets in a tab and check if all the radiobuttons work properly
    """
    for widget in widgets:
        try:
            radio_button_x_axis = widget.tk_radio_button_1
            radio_button_STD_axis = widget.tk_radio_button_2
        # The first part of the loop makes sure the buttons work individually
            
            # This part checks for the x-axis radiobuttons
            radio_button_x_axis.invoke()
            assert radio_button_x_axis['value'] == gui.x_axis_label_widget.get()
            gui.x_axis_label_widget.set(None)

            # This part checks for the std-axis radiobuttons
            radio_button_STD_axis.invoke()
            assert radio_button_STD_axis['value'] == gui.STD_axis_label_widget.get()
            gui.x_axis_label_widget.set(None)
            
                

        # The second part of the loop checks if correct ValueErrors get raised if both selected axis are the same
            # This part necessitated that there are no message boxes
            radio_button_x_axis.invoke()
            radio_button_STD_axis.invoke()
            with raises(ValueError):
                gui.update_plot()
        except:
            assert widget.can_be_axis == False


def test_source_widgets() -> None:
    widgets = gui.source_tab.widget_collection.widgets
    radio_button_invoke(widgets)


def test_scatterer_widgets() -> None:
    for tab in gui.scatterer_tab.type_widget['values']:
        gui.scatterer_tab.type_widget.set(tab)
        gui.scatterer_tab.on_type_change()
        widgets = gui.scatterer_tab.widget_collection.widgets
        radio_button_invoke(widgets)


def test_detector_widgets() -> None:
    for tab in gui.detector_tab.type_widget['values']:
        gui.detector_tab.type_widget.set(tab)
        gui.detector_tab.on_type_change()
        widgets = gui.detector_tab.widget_collection.widgets
        radio_button_invoke(widgets)


