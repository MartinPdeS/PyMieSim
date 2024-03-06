import tkinter as tk
from tkinter import ttk, filedialog
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# Assuming PyMieSim and its dependencies are installed
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment import Setup
from PyMieSim import measure


class PyMieSimGUI:
    """
    A graphical user interface for computing and visualizing the B1 scattering coefficient
    for cylindrical scatterers using PyMieSim.

    Attributes:
        master: The main tkinter window.
    """

    def __init__(self, master):
        """
        Initializes the GUI, setting up widgets and internal state.
        """
        self.master = master
        self.master.title("PyMieSim Graphic Interface")

        self.vars = self.setup_variables()
        self.vars_dict = {label: var for label, var in self.vars}
        self.setup_widgets()

    def setup_variables(self):
        """
        Setup and return a list of tuples containing label names and associated tkinter StringVars.
        """
        return [
            ('Wavelength (nm)', tk.StringVar(value='400')),
            ('Polarization Value (radians)', tk.StringVar(value='0')),
            ('Optical Power (mW)', tk.StringVar(value='1')),
            ('Numerical Aperture (NA)', tk.StringVar(value='0.2')),
            ('Diameter Start (nm)', tk.StringVar(value='100')),
            ('Diameter End (nm)', tk.StringVar(value='10000')),
            ('Diameter Points', tk.StringVar(value='800')),
            ('Index', tk.StringVar(value='1.4')),
            ('Medium Refractive Index', tk.StringVar(value='1')),
        ]

    def setup_widgets(self):
        """
        Sets up the input fields, buttons, and plot frame within the GUI.
        """
        self.inputs_frame = ttk.Frame(self.master)
        self.inputs_frame.pack(side=tk.LEFT, fill=tk.Y)

        self.plot_frame = tk.Frame(self.master)
        self.plot_frame.pack(fill=tk.BOTH, expand=1)

        for i, (label, var) in enumerate(self.vars):
            ttk.Label(self.inputs_frame, text=label).grid(column=0, row=i, sticky=tk.W)
            ttk.Entry(self.inputs_frame, textvariable=var).grid(column=1, row=i)

        ttk.Button(self.inputs_frame, text="Calculate", command=self.update_plot).grid(column=0, row=len(self.vars), columnspan=2, pady=10)
        ttk.Button(self.inputs_frame, text="Save as CSV", command=self.save_data_as_csv).grid(column=0, row=len(self.vars) + 1, columnspan=2, pady=10)

    def get_data_from_gui(self):
        """
        Retrieves input values from the GUI, sets up PyMieSim components, and computes the B1 scattering data.
        """
        self.gui_wavelength = float(self.vars_dict['Wavelength (nm)'].get()) * 1e-9
        self.gui_polarization_value = float(self.vars_dict['Polarization Value (radians)'].get())
        self.gui_optical_power = float(self.vars_dict['Optical Power (mW)'].get()) * 1e-3
        self.gui_NA = float(self.vars_dict['Numerical Aperture (NA)'].get())
        self.gui_diameter_start = float(self.vars_dict['Diameter Start (nm)'].get()) * 1e-9
        self.gui_diameter_end = float(self.vars_dict['Diameter End (nm)'].get()) * 1e-9
        self.gui_diameter_points = int(self.vars_dict['Diameter Points'].get())
        self.gui_index = float(self.vars_dict['Index'].get())
        self.gui_n_medium = float(self.vars_dict['Medium Refractive Index'].get())

    def get_data_from_PyMieSim(self) -> None:
        self.source_set = Gaussian(
            wavelength=self.gui_wavelength,
            polarization_value=self.gui_polarization_value,
            polarization_type='linear',
            optical_power=self.gui_optical_power,
            NA=self.gui_NA
        )

        self.scatterer_set = Cylinder(
            diameter=np.linspace(self.gui_diameter_start, self.gui_diameter_end, self.gui_diameter_points),
            index=self.gui_index,
            n_medium=self.gui_n_medium,
            source_set=self.source_set
        )

        self.experiment = Setup(
            scatterer_set=self.scatterer_set,
            source_set=self.source_set
        )

        self.data = self.experiment.get(measure.Qsca)

        return self.data

    def save_data_as_csv(self):
        """
        Triggered by the "Save as CSV" button. Opens a file dialog to save the computed data as a CSV file.
        """
        if hasattr(self, 'data'):
            filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
            if filepath:
                # Assuming self.data is a pandas DataFrame or can be converted to one
                # Adjust this part based on the actual data structure of self.data
                self.data.to_csv(filepath, index=False)
                print(f"Data saved to {filepath}")
        else:
            print("No data to save. Please calculate first.")

    def update_plot(self):
        """
        Triggered by the "Calculate" button. It computes the coefficient and updates the plot.
        """
        self.get_data_from_gui()

        self.get_data_from_PyMieSim()

        figure = self.data.plot(x=self.scatterer_set.diameter)
        figure.unit_size = (7, 4)
        figure._render_()
        figure = figure._mpl_figure

        # plot.plot(self.scatterer_set.diameter.values * 1e9, self.data, label='B1 Scattering Coefficient')

        # Embed the plot
        for widget in self.plot_frame.winfo_children():
            widget.destroy()

        canvas = FigureCanvasTkAgg(figure, master=self.plot_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


def main():
    root = tk.Tk()
    PyMieSimGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()

# -
