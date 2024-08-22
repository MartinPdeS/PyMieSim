#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tkinter.ttk import Frame
from tkinter import messagebox, filedialog
import tkinter
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from typing import NoReturn, Dict

from PyMieSim.gui.singleton import datashelf
from PyMieSim.experiment import Setup
from PyMieSim.gui.widget_collection import WidgetCollection


class ControlTab:
    """
    This class is defines and handles the control buttons. It define the following:
        - calculate_button: used to produce a graph
        - save_button: used to save the data from the graph as a csv fils
        - export_button: used to export the graph
        - reset_std_button: used to undo the std axis selection on the gui

    Other attributes:
    -frame (ttk.Frame)
    """

    def __init__(self, frame: Frame) -> None:
        self.frame = frame
        self.setup_widgets()

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

    def setup_widgets(self):
        """"
        Creates the control buttons, and binds them to their specific command
        """

        # Dictonary with the following format: {"button_nam": button_command}
        self.button_config = {
            "calculate_button": self.calculate_plot,
            "save_button": self.save_data_as_csv,
            "export_button": self.export_plot,
            "reset_std_button": self.reset_STDaxis_selection
        }

        self.widget_collection = WidgetCollection(frame=self.frame)

        self.widget_collection.setup_control_widget(config_dict=self.button_config)

        for widget in self.widget_collection.widgets:
            setattr(self, widget.component_label, widget)

    def save_data_as_csv(self) -> NoReturn:
        """
        Triggered by the "Save as CSV" button. Opens a file dialog to save the computed data as a CSV file.
        """

        if hasattr(datashelf, 'data'):
            self.filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
            if self.filepath:

                # Assuming datashelf.data is a pandas DataFrame or can be converted to one
                numpy.savetxt(self.filepath, datashelf.data.y.values.squeeze(), delimiter=",")
                print(f"Data saved to {self.filepath}")
        else:
            print("No data to save. Please calculate first.")

    def export_plot(self) -> NoReturn:
        """
        Opens a file dialog for the user to choose where to save the current plot,
        then saves the plot to the specified location.
        """
        # Ensure there's a plot to save
        if hasattr(datashelf, 'figure'):
            # Open file dialog to choose file name and type
            filetypes = [
                ('PNG files', '*.png'),
                ('JPEG files', '*.jpg;*.jpeg'),
                ('PDF files', '*.pdf'),
                ('SVG files', '*.svg'),
                ('All files', '*.*')
            ]

            self.filepath = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=filetypes,
                title="Save plot as..."
            )

            # If a file was selected (i.e., dialog not cancelled)
            if self.filepath:
                # Save the figure using matplotlib's savefig
                datashelf.figure.savefig(self.filepath)
                messagebox.showinfo("Export Successful", f"Plot successfully saved to {self.filepath}")
        else:
            messagebox.showwarning("Export Failed", "No plot available to export.")

    def reset_STDaxis_selection(self):
        """
        Allows the user to unselect the std-axis radiobuttons.
        """
        datashelf.STD_axis_label_widget.set(None)

    def setup_experiment(self) -> NoReturn:
        """
        Compute the B1 scattering data using either a single diameter or a range of diameters.
        """
        datashelf.scatterer_tab.setup_component()
        datashelf.source_tab.setup_component()
        datashelf.detector_tab.setup_component()

        self.experiment = Setup(
            scatterer=datashelf.scatterer_tab.component,
            source=datashelf.source_tab.component,
            detector=datashelf.detector_tab.component
        )

    def validate_axis_choice(self):
        if self.x_axis == self.std_axis:
            return "Warning: x-axis cannot be equal to STD-axis."

        if self.y_axis_selection != "coupling" and self.std_axis in self.detector_tab.component_dict.keys():
            return "Warning: STD-axis cannot be associated to detector if y-axis is not coupling."

        if self.y_axis_selection != "coupling" and self.x_axis in self.detector_tab.component_dict.keys():
            return "Warning: x-axis cannot be associated to detector if y-axis is not coupling."

        return True

    def calculate_plot(self) -> NoReturn:
        # Closing all previous plots
        plt.close('all')

        # Defining the axis'
        self.x_axis = datashelf.x_axis_label_widget.get()
        self.std_axis = datashelf.STD_axis_label_widget.get()
        self.y_axis_selection = datashelf.y_axis_selection.get()

        # Checking if axis selection is valid
        validation = self.validate_axis_choice()
        if validation is not True:
            self.messagebox = messagebox.showerror(title="error", message=validation, parent=self.frame)
            raise ValueError(validation)

        # Setting up the data and the components
        y_axis = datashelf.measure_map[self.y_axis_selection]

        self.setup_experiment()

        datashelf.data = self.experiment.get(y_axis)

        self.x_axis_component = self.axis_mapping[self.x_axis]

        self.STD_axis_component = None if self.std_axis == "None" else self.axis_mapping[self.std_axis]

        try:
            self.generate_figure()

        except ValueError as e:
            messagebox.showerror("Input Error", str(e))

    def generate_figure(self) -> NoReturn:
        """
        Generates and displays the simulation results as a plot in a new window.
        """
        if hasattr(self, 'new_window'):
            self.new_window.destroy()

        # Creates a tk window for the plot
        self.new_window = tkinter.Toplevel(self.frame)
        self.new_window.title("Plot Window")

        # Renders the figure
        figure = datashelf.data.plot(x=self.x_axis_component, std=self.STD_axis_component)
        figure.unit_size = (9, 4)
        figure._render_()
        datashelf.figure = figure._mpl_figure

        # Creates the canvas
        canvas = FigureCanvasTkAgg(datashelf.figure, master=self.new_window)
        canvas.get_tk_widget().pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(canvas, self.new_window)
        canvas._tkcanvas.pack(side=tkinter.TOP, fill=tkinter.BOTH, expand=True)

        # Draws the figure
        canvas.draw()
        self.toolbar.update()

# -
