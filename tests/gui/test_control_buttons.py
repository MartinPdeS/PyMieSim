#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import random
import tkinter
from PyMieSim.gui.main_window import PyMieSimGUI
from PyMieSim.experiment.measure import __sphere__
from unittest.mock import patch, MagicMock
import DataVisual

def random_list(min, max):
    lenght = 3
    return [random.randint(min, max) for i in range(lenght)]


def reset_std_button(widget, gui):
    '''
    This is the function that will assert whether the reset_std_button works properly
    '''
    gui.reset_std_button.invoke()
    assert gui.STD_axis_label_widget.get() == 'None', f"reset button did not worlkfor {widget.tk_radio_button_2['value']} widget"


@patch('DataVisual.multi_array.Array.plot')
@patch('tkinter.messagebox.showerror')
def calculate_button(mock_messagebox, mock_draw, gui, y_axis_index, source_widgets, scatterer_widgets, detector_widgets):
    """
    This function ensures that the calculate button generates a graph for
    all possible x-axis values, along with some combination of random std-axis values.
    It does not test for all possible combinations of x and std axis due to computational time constraints.
    """

    # Axis choices from the detector tab only makes sense if the y_axis is coupling
    possible_widgets = [*source_widgets, *scatterer_widgets]
    if y_axis_index == 21:
        possible_widgets.extend(detector_widgets)


    for x_widget in possible_widgets:
        if x_widget.tk_radio_button_1['value'] == 'mode_number':
            continue
        mock_draw.reset_mock()

        x_widget.tk_radio_button_1.invoke()
        gui.calculate_button.invoke()

        assert mock_draw.call_count == 1, f"calculate_button with x-axis selection '{x_widget.tk_radio_button_1['value']}' did not call the draw"

        for position in random_list(0,len(possible_widgets)-1): # This will make a test for 3 random std_axis selection (doing them all takes too much computational time)
            mock_draw.reset_mock()

            std_widget = possible_widgets[position]
            std_widget.tk_radio_button_2.invoke()
            gui.calculate_button.invoke()

            try:
                assert mock_draw.call_count == 1
            except:
                assert x_widget.tk_radio_button_1['value'] == std_widget.tk_radio_button_2['value'], f"calculate_button with x-axis selection '{x_widget.tk_radio_button_1['value']}' and std-axis selection '{std_widget.tk_radio_button_2['value']}' did not call the draw as intended"

            reset_std_button(widget=std_widget, gui=gui)

@pytest.mark.parametrize('y_axis_index', [index for index in range(len(__sphere__))], ids=__sphere__) 
def test_in_all_combination_of_widgets(y_axis_index):
    """
    This function is meant to cycle through all combinations of tabs and generate the corresponding battery of 
    widgets with which the tests should be run. It will then execute the tests for all tab combinations.
    """
    root = tkinter.Tk()
    root.geometry("750x600")
    gui = PyMieSimGUI(root)

    y_axis_widget = gui.axis_tab.widget_collection.widgets[0]
    y_axis_widget.tk_widget.current(y_axis_index)
    
    measure_input = gui.axis_tab.get_inputs()[0]
    if measure_input not in ['Qsca', 'Csca', 'Qext', 'Cext', 'Qabs', 'Cabs', 'coupling'] :
        root.destroy()
        return

    # The following nested for loops will create all possible widget combinations

    source_widgets = gui.source_tab.widget_collection.widgets[0:2] # Only the first two widgets can be choosen as axis

    for scatterer_str in gui.scatterer_tab.type_widget['values']:
        gui.scatterer_tab.type_widget.set(scatterer_str)
        gui.scatterer_tab.on_type_change()
        scatterer_widgets = gui.scatterer_tab.widget_collection.widgets

        for detector_str in gui.detector_tab.type_widget['values']:
            gui.detector_tab.type_widget.set(detector_str)
            gui.detector_tab.on_type_change()
            all_detector_widgets = gui.detector_tab.widget_collection.widgets
            detector_widgets =[widget for widget in all_detector_widgets if widget.component_label != "mean_coupling"]

            calculate_button(
                gui=gui,
                y_axis_index=y_axis_index,
                source_widgets=source_widgets,
                scatterer_widgets=scatterer_widgets,
                detector_widgets=detector_widgets
            )

    root.destroy()

@patch('tkinter.filedialog.asksaveasfilename')
@patch('tkinter.messagebox.showinfo')
def test_export_plot_button(mock_messagebox, mock_filepath):
    """
    This function tests whether the export_plot_button function is working as intended.
    """
    
    # setting up the environment
    root = tkinter.Tk()
    root.geometry("750x600")
    gui = PyMieSimGUI(root)

    # mocking the necessary variables
    gui.figure = MagicMock()
    gui.filepath = MagicMock()

    # invoking the button
    gui.export_button.invoke()

    # the assertion
    assert gui.figure.savefig.call_count == 1

    root.destroy()


@patch('numpy.savetxt')
@patch('tkinter.filedialog.asksaveasfilename')
def test_save_as_csv_button(mock_filepath, mock_save):
    """
    This test takes a battery of widgets and checks if the calculate button calls 
    the draw method (i.e., creates a graph) for all selections of x-axis buttons.
    """

    # setting up the environment
    root = tkinter.Tk()
    root.geometry("750x600")
    gui = PyMieSimGUI(root)

    # mocking the necessary variables
    gui.data = MagicMock()
    gui.filepath = MagicMock()

    # invoking the button
    gui.save_button.invoke()

    # the assertion
    assert mock_save.call_count == 1
    
    root.destroy()
    
