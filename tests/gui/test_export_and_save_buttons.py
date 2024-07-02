#!/usr/bin/env python
# -*- coding: utf-8 -*-

import tkinter

from PyMieSim.gui.main_window import PyMieSimGUI
from PyMieSim.gui.singleton import datashelf

from unittest.mock import patch, MagicMock


@patch('tkinter.filedialog.asksaveasfilename')
@patch('tkinter.messagebox.showinfo')
def test_export_plot_button(mock_messagebox, mock_filepath):
    """
    This function tests whether the export_plot_button function is working as intended.
    """

    # setting up the environment
    root = tkinter.Tk()
    root.geometry("750x600")
    PyMieSimGUI(root)

    # mocking the necessary variables
    datashelf.figure = MagicMock()
    datashelf.control_tab.filepath = MagicMock()

    # invoking the button
    datashelf.control_tab.export_button.invoke()

    # the assertion
    assert datashelf.figure.savefig.call_count == 1

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
    PyMieSimGUI(root)

    # mocking the necessary variables
    datashelf.data = MagicMock()
    datashelf.control_tab.filepath = MagicMock()

    # invoking the button
    datashelf.control_tab.save_button.invoke()

    # the assertion
    assert mock_save.call_count == 1

    root.destroy()

# -
