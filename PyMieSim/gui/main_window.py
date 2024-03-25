#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
# Assuming PyMieSim and its dependencies are installed
from PyMieSim.experiment import Setup

from PyMieSim.gui import SourceTab, ScattererTab, DetectorTab, AxisTab


class PyMieSimGUI:
    """
    Graphical User Interface for computing and visualizing the B1 scattering coefficient
    for cylindrical scatterers using PyMieSim.

    Attributes:
        master (tk.Tk): The main tkinter window.
    """

    def __init__(self, master: tk.Tk):
        """
        Initializes the GUI, setting up variables, plot frame, notebook, and controls.

        Parameters:
            master (tk.Tk): The root window of the application.
        """
        self.master = master
        self.master.protocol("WM_DELETE_WINDOW", self.on_close)
        self.master.title("PyMieSim Graphic Interface")

        self.customize_notebook_style()
        self.setup_plot_frame()
        self.setup_notebook()
        self.setup_controls()

    def on_close(self):
        """Handles the GUI close event."""
        plt.close('all')  # Close all matplotlib figures
        self.master.destroy()  # Close the Tkinter window

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

        # Style for larger buttons
        style.configure(
            "Large.TButton",
            font=('Helvetica', 18),    # Font size for buttons
            padding=[20, 15]           # Increased padding for top and bottom to make the button taller
        )

    def setup_notebook(self) -> NoReturn:
        """
        Sets up the notebook widget with tabs for Source, Scatterer, and Detector configurations.
        """
        self.notebook = ttk.Notebook(self.master)
        self.notebook.pack(fill=tk.BOTH, expand=True, side=tk.TOP)

        # Create tab instances
        self.source_tab = SourceTab(notebook=self.notebook, label='Source')
        self.scatterer_tab = ScattererTab(self.notebook, 'Scatterer', source_tab=self.source_tab, main_window=self)
        self.detector_tab = DetectorTab(self.notebook, 'Detector')
        self.axis_tab = AxisTab(self.notebook, 'Axis Configuration', other_tabs=[self.source_tab, self.scatterer_tab, self.detector_tab])

    def export_plot(self) -> NoReturn:
        """
        Opens a file dialog for the user to choose where to save the current plot,
        then saves the plot to the specified location.
        """
        # Ensure there's a plot to save
        if hasattr(self, 'figure'):
            # Open file dialog to choose file name and type
            filetypes = [
                ('PNG files', '*.png'),
                ('JPEG files', '*.jpg;*.jpeg'),
                ('PDF files', '*.pdf'),
                ('SVG files', '*.svg'),
                ('All files', '*.*')
            ]

            filepath = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=filetypes,
                title="Save plot as..."
            )

            # If a file was selected (i.e., dialog not cancelled)
            if filepath:
                # Save the figure using matplotlib's savefig
                self.figure.savefig(filepath)
                messagebox.showinfo("Export Successful", f"Plot successfully saved to {filepath}")
        else:
            messagebox.showwarning("Export Failed", "No plot available to export.")

    def setup_controls(self) -> NoReturn:
        """
        Sets up control buttons for calculating results and saving data.
        """
        self.controls_frame = ttk.Frame(self.master)
        self.controls_frame.pack(fill=tk.X, side=tk.TOP)

        ttk.Button(
            self.controls_frame,
            text="Calculate",
            style="Large.TButton",  # Apply the custom style for larger buttons
            command=self.update_plot
        ).pack(side=tk.LEFT, fill=tk.BOTH)

        ttk.Button(
            self.controls_frame,
            text="Save as CSV",
            style="Large.TButton",  # Apply the custom style for larger buttons
            command=self.save_data_as_csv
        ).pack(side=tk.LEFT, fill=tk.BOTH)

        # Add Export Plot button
        ttk.Button(
            self.controls_frame,
            text="Export Plot",
            style="Large.TButton",
            command=self.export_plot  # Method to be implemented
        ).pack(side=tk.LEFT, fill=tk.BOTH)

    def setup_plot_frame(self) -> NoReturn:
        """
        Sets up the frame for displaying plots.
        """
        self.plot_frame = tk.Frame(self.master)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

    def setup_PyMieSim(self) -> NoReturn:
        """
        Compute the B1 scattering data using either a single diameter or a range of diameters.
        """
        self.source_tab.setup_component()
        self.scatterer_tab.setup_component()
        self.detector_tab.setup_component()

        self.experiment = Setup(
            scatterer_set=self.scatterer_tab.component,
            source_set=self.source_tab.component,
            detector_set=self.detector_tab.component
        )

    def save_data_as_csv(self) -> NoReturn:
        """
        Triggered by the "Save as CSV" button. Opens a file dialog to save the computed data as a CSV file.
        """
        if hasattr(self, 'data'):
            filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
            if filepath:
                # Assuming self.data is a pandas DataFrame or can be converted to one
                np.savetxt(filepath, self.data.y.values.squeeze(), delimiter=",")
                print(f"Data saved to {filepath}")
        else:
            print("No data to save. Please calculate first.")

    def generate_figure(self, x_axis: object) -> NoReturn:
        figure = self.data.plot(x=x_axis)
        figure.unit_size = (9, 4)
        figure._render_()
        self.figure = figure._mpl_figure

        # Embed the plot
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(self.figure, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Add the navigation toolbar
        self.toolbar = NavigationToolbar2Tk(canvas, self.plot_frame)
        self.toolbar.update()
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def update_plot(self) -> NoReturn:
        plt.close('all')
        x_axis, y_axis = self.axis_tab.get_inputs()

        self.y_axis = self.axis_tab.measure_map[y_axis]

        self.setup_PyMieSim()

        self.data = self.experiment.get(self.y_axis)

        self.x_axis = self.axis_tab.axis_mapping['wavelength']

        try:
            self.generate_figure(self.axis_tab.x_axis)

        except ValueError as e:
            messagebox.showerror("Input Error", str(e))


# -
