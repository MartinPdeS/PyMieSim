#!/usr/bin/env python
# -*- coding: utf-8 -*-


from typing import NoReturn
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

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
        self.link_radio_button = "link"
        self.x_axis_label_widget = tk.StringVar(value='phi_offset')
        self.STD_axis_label_widget = tk.StringVar(value=None)
        self.STD_axis_label_widget.set(None)
        self.customize_notebook_style()
        self.setup_notebook()
        self.setup_controls()

    def on_close(self) -> NoReturn:
        """
        Handles the GUI close event.
        """
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

        style.configure(
            "Large.TButton",
            font=('Helvetica', 18),
            padding=[20, 20]
        )

    def setup_notebook(self) -> NoReturn:
        """
        Sets up the notebook widget with tabs for Source, Scatterer, and Detector configurations.
        """
        self.notebook = ttk.Notebook(self.master)
        self.notebook.grid(row=0, column=0, sticky="ewns")

        self.notebook_2 = ttk.Notebook(self.master)
        self.notebook_2.grid(row=2, column=0, sticky="ewns")

        # Create tab instances
        self.source_tab = SourceTab(
            x_axis=self.x_axis_label_widget,
            STD_axis=self.STD_axis_label_widget,
            notebook=self.notebook,
            label='Source'
        )

        self.scatterer_tab = ScattererTab(
            x_axis=self.x_axis_label_widget,
            STD_axis=self.STD_axis_label_widget,
            notebook=self.notebook,
            label='Scatterer',
            source_tab=self.source_tab
        )

        self.detector_tab = DetectorTab(
            x_axis=self.x_axis_label_widget,
            STD_axis=self.STD_axis_label_widget,
            notebook=self.notebook,
            label='Detector'
        )

        self.axis_tab = AxisTab(
            notebook=self.notebook_2,
            label='Axis Configuration',
            other_tabs=[self.source_tab, self.scatterer_tab, self.detector_tab]
        )

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

            self.filepath = filedialog.asksaveasfilename(
                defaultextension=".png",
                filetypes=filetypes,
                title="Save plot as..."
            )

            # If a file was selected (i.e., dialog not cancelled)
            if self.filepath:
                # Save the figure using matplotlib's savefig
                self.figure.savefig(self.filepath)
                messagebox.showinfo("Export Successful", f"Plot successfully saved to {self.filepath}")
        else:
            messagebox.showwarning("Export Failed", "No plot available to export.")

    def setup_controls(self) -> NoReturn:
        """
        Sets up control buttons for calculating results and saving data.
        """
        self.controls_frame = ttk.Frame(self.master)
        self.controls_frame.grid(row=1, column=0, sticky="ew")

        self.calculate_button = ttk.Button(
            self.controls_frame,
            text="Calculate",
            style="Large.TButton",
            command=self.update_plot
        )
        self.calculate_button.grid(row=0, column=0, sticky="ew")

        self.save_button = ttk.Button(
            self.controls_frame,
            text="Save as CSV",
            style="Large.TButton",
            command=self.save_data_as_csv
        )
        self.save_button.grid(row=0, column=1, sticky="ew")

        self.export_button = ttk.Button(
            self.controls_frame,
            text="Export Plot",
            style="Large.TButton",
            command=self.export_plot
        )
        self.export_button.grid(row=0, column=2, sticky="ew")

        self.reset_std_button = ttk.Button(
            self.controls_frame,
            text="Reset STD-axis",
            style="Large.TButton",
            command=self.reset_STDaxis_selection
        )
        self.reset_std_button.grid(row=0, column=3, sticky="ew")

    def reset_STDaxis_selection(self):
        """
        Allows the user to unselect the std-axis radiobuttons.
        """
        self.STD_axis_label_widget.set(None)

    def setup_PyMieSim(self) -> NoReturn:
        """
        Compute the B1 scattering data using either a single diameter or a range of diameters.
        """
        self.scatterer_tab.setup_component()
        self.source_tab.setup_component()
        self.detector_tab.setup_component()

        self.experiment = Setup(
            scatterer=self.scatterer_tab.component,
            source=self.source_tab.component,
            detector=self.detector_tab.component
        )

    def save_data_as_csv(self) -> NoReturn:
        """
        Triggered by the "Save as CSV" button. Opens a file dialog to save the computed data as a CSV file.
        """

        if hasattr(self, 'data'):
            self.filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
            if self.filepath:
                # Assuming self.data is a pandas DataFrame or can be converted to one
                np.savetxt(self.filepath, self.data.y.values.squeeze(), delimiter=",")
                print(f"Data saved to {self.filepath}")
        else:
            print("No data to save. Please calculate first.")

    def generate_figure(self):
        """
        Generates and displays the simulation results as a plot in a new window.
        """
        if hasattr(self, 'new_window'):
            self.new_window.destroy()

        self.new_window = tk.Toplevel(self.master)
        self.new_window.title("Plot Window")

        figure = self.data.plot(x=self.x_axis_component, std=self.STD_axis_component)
        figure.unit_size = (9, 4)
        figure._render_()
        self.figure = figure._mpl_figure

        canvas = FigureCanvasTkAgg(self.figure, master=self.new_window)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.toolbar = NavigationToolbar2Tk(canvas, self.new_window)
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        canvas.draw()
        self.toolbar.update()

    def update_plot(self) -> NoReturn:
        plt.close('all')

        x_axis, y_axis, std_axis = self.x_axis_label_widget.get(), self.axis_tab.get_inputs()[0], self.STD_axis_label_widget.get()

        if x_axis == std_axis:
            self.messagebox1 = tk.messagebox.showerror(title="error", message="X-axis cannot be equal to STD-axis.", parent=self.master)
            raise ValueError("Warning: x-axis cannot be equal to STD-axis.")

        if y_axis != "coupling" and std_axis in self.detector_tab.component_dict.keys():
            self.messagebox2 = tk.messagebox.showerror(title="error", message="STD-axis cannot be associated to detector if y-axis is not coupling.", parent=self.master)
            raise ValueError("Warning: STD-axis cannot be associated to detector if y-axis is not coupling.")

        if y_axis != "coupling" and x_axis in self.detector_tab.component_dict.keys():
            self.messagebox3 = tk.messagebox.showerror(title="error", message="x-axis cannot be associated to detector if y-axis is not coupling.", parent=self.master)
            raise ValueError("Warning: x-axis cannot be associated to detector if y-axis is not coupling.")

        self.y_axis = self.axis_tab.measure_map[y_axis]

        self.setup_PyMieSim()

        self.data = self.experiment.get(self.y_axis)

        self.x_axis_component = self.axis_tab.axis_mapping[x_axis]

        self.STD_axis_component = None if std_axis == "None" else self.axis_tab.axis_mapping[std_axis]

        try:
            self.generate_figure()

        except ValueError as e:
            messagebox.showerror("Input Error", str(e))

# -
