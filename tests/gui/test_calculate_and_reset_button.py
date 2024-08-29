#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest
import tkinter
import itertools

from PyMieSim.gui.main_window import PyMieSimGUI
from PyMieSim.gui.singleton import datashelf
from unittest.mock import patch


measures = ['Qsca', 'Csca', 'Qabs', 'coupling']  # A selection of measures to test


@patch('tkinter.messagebox.showerror')
@patch('PyMieSim.gui.control_tab.ControlTab.generate_figure')
def calculate_and_reset_button(mock_plot, mock_message_box, possible_widgets):
    """
    This function ensures that the calculate button generates a graph for
    all possible x-axis vues, along with some combination of random std-axis values.
    It does not test for all possible combinations of x and std axis due to computational time constraints.
    """
    for (x_widget, std_widget) in itertools.combinations(possible_widgets, 2):
        mock_plot.reset_mock()

        if x_widget.tk_radio_button_1['value'] == 'mode_number':
            continue

        x_widget.tk_radio_button_1.invoke()
        datashelf.control_tab.calculate_button.invoke()

        assert mock_plot.call_count == 1, f"calculate_button with x-axis selection '{x_widget.tk_radio_button_1['value']}' did not call the draw"

        std_widget.tk_radio_button_2.invoke()
        datashelf.control_tab.calculate_button.invoke()

        datashelf.control_tab.reset_std_button.invoke()
        assert datashelf.STD_axis_label_widget.get() == 'None', f"reset button did not worlkfor {std_widget.tk_radio_button_2['value']} widget"


@pytest.mark.parametrize('measure', measures)
@pytest.mark.parametrize('detector_tab', ['Photodiode', 'CoherentMode'])
@pytest.mark.parametrize('scatterer_tab', ['Sphere', 'Cylinder', 'CoreShell'])
def test_in_all_combination_of_widgets(scatterer_tab, detector_tab, measure):
    """
    This function is meant to cycle through all combinations of tabs and generate the corresponding battery of
    widgets with which the tests should be run. It will then execute the tests for all tab combinations.
    """

    # Setting up the GUI
    root = tkinter.Tk()
    root.geometry("750x600")
    PyMieSimGUI(root)

    # Choose the correct measure (i.e. y-axis)
    datashelf.axis_tab.widget_collection.widgets[0].tk_widget.set(measure)

    # Set up the tabs
    datashelf.scatterer_tab.combobox_widget_collection.combobox_widget.tk_widget.set(scatterer_tab)
    datashelf.scatterer_tab.on_type_change()

    datashelf.detector_tab.combobox_widget_collection.combobox_widget.tk_widget.set(detector_tab)
    datashelf.detector_tab.on_type_change()

    # The widgets collections
    source_widgets = datashelf.source_tab.widget_collection
    scatterer_widgets = datashelf.scatterer_tab.widget_collection
    detector_widgets = datashelf.detector_tab.widget_collection

    possible_widgets = [
        source_widgets['wavelength'],
        scatterer_widgets['medium_index'],
        detector_widgets['NA'],
        detector_widgets['gamma_offset'],
        detector_widgets['polarization_filter']
    ]

    # Run the test
    calculate_and_reset_button(possible_widgets=possible_widgets)

    root.destroy()

if __name__ == "__main__":
    pytest.main([__file__])

# -
