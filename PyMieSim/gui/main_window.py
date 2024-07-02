#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import NoReturn, Dict
import tkinter
from tkinter import ttk
import matplotlib.pyplot as plt

from PyMieSim.gui.singleton import datashelf
from PyMieSim.gui import SourceTab, ScattererTab, DetectorTab, AxisTab
from PyMieSim.gui.control_tab import ControlTab


class PyMieSimGUI:
    """
    Graphical User Interface for computing and visualizing the B1 scattering coefficient
    for cylindrical scatterers using PyMieSim.

    Attributes:
        master (tkinter.Tk): The main tkinter window.
    """

    def __init__(self, master: tkinter.Tk):
        """
        Initializes the GUI, setting up variables, plot frame, notebook, and controls.

        Parameters:
            master (tk.Tk): The root window of the application.
        """
        self.master = master
        self.master.protocol("WM_DELETE_WINDOW", self.on_close)
        self.master.title("PyMieSim Graphic Interface")

        datashelf.x_axis_label_widget = tkinter.StringVar(value='phi_offset')
        datashelf.STD_axis_label_widget = tkinter.StringVar(value=None)
        datashelf.STD_axis_label_widget.set(None)

        self.link_radio_button = "link"
        self.customize_notebook_style()
        self.setup_notebook()

    def on_close(self) -> NoReturn:
        """
        Handles the GUI close event.
        """
        plt.close('all')  # Close all matplotlib figures
        self.master.destroy()  # Close the Tkinter window

    @property
    def axis_mapping(self) -> Dict[str, str]:
        """
        Combines mappings from all other tabs to provide a comprehensive dictionary of available axis options.

        Returns:
            Dict[str, str]: A dictionary mapping UI labels to internal scatterer parameter names.
        """
        _axis_mapping = {}
        for tab in [datashelf.source_tab, datashelf.detector_tab, datashelf.scatterer_tab]:
            _axis_mapping.update(tab.component.mapping)

        return _axis_mapping

    # The following section of the class will setup the notebooks and their content
    def setup_notebook(self) -> NoReturn:
        """
        Sets up the notebook widget with tabs for Source, Scatterer, and Detector configurations.
        """
        self.notebook = ttk.Notebook(self.master)
        self.notebook.grid(row=0, column=0, sticky="ewns")

        self.notebook_2 = ttk.Notebook(self.master)
        self.notebook_2.grid(row=2, column=0, sticky="ewns")

        self.controls_frame = ttk.Frame(self.master)
        self.controls_frame.grid(row=1, column=0, sticky="ew")

        # Create tab instances, in datashelf, for easy access by other classes
        datashelf.axis_tab = AxisTab(
            notebook=self.notebook_2,
            label='Axis Configuration'
        )

        datashelf.control_tab = ControlTab(
            frame=self.controls_frame,
        )

        datashelf.source_tab = SourceTab(
            notebook=self.notebook,
            label='Source'
        )

        datashelf.scatterer_tab = ScattererTab(
            notebook=self.notebook,
            label='Scatterer',
            axis_tab=datashelf.axis_tab
        )

        datashelf.detector_tab = DetectorTab(
            notebook=self.notebook,
            label='Detector'
        )

    def customize_notebook_style(self) -> NoReturn:
        """
        Customizes the ttk Notebook style for a unique appearance of tabs, making them larger.
        """
        style = ttk.Style()
        style.configure("TNotebook", background="#f0f0f0")

        style.configure(
            "TNotebook.Tab",
            background="#d0d0d0",
            padding=[10, 20, 10, 20],  # Increase padding for larger tabs
            font=('Helvetica', 12)     # Larger font for tabs
        )

        style.map(
            "TNotebook.Tab",
            background=[("selected", "#a0a0a0")],
            expand=[("selected", [1, 1, 1, 0])]
        )

        style.configure(
            "Large.TButton",
            font=('Helvetica', 18),
            padding=[20, 20]
        )

# -
