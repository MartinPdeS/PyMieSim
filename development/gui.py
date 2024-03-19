import tkinter as tk
from tkinter import ttk, filedialog, messagebox
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
        self.master = master
        self.master.title("PyMieSim Graphic Interface")
        self.setup_variables()
        self.setup_plot_frame()
        self.setup_notebook()

    def setup_variables(self):
        self.source_vars = {
            'Wavelength (nm)': tk.StringVar(value='400'),
            'Polarization Value (radians)': tk.StringVar(value='0'),
            'Optical Power (mW)': tk.StringVar(value='1'),
            'Numerical Aperture (NA)': tk.StringVar(value='0.2'),
        }

        self.scatterer_vars = {
            'Diameter (nm)': tk.StringVar(value='100'),  # Accepts both single value and range
            'Index': tk.StringVar(value='1.4'),
            'Medium Refractive Index': tk.StringVar(value='1'),
        }

    def setup_notebook(self):
        self.notebook = ttk.Notebook(self.master)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        self.source_tab = ttk.Frame(self.notebook)
        self.scatterer_tab = ttk.Frame(self.notebook)

        self.notebook.add(self.source_tab, text='Source')
        self.notebook.add(self.scatterer_tab, text='Scatterer')

        self.setup_source_tab()
        self.setup_scatterer_tab()

    def setup_source_tab(self):
        row = 0
        for label, var in self.source_vars.items():
            ttk.Label(self.source_tab, text=label).grid(column=0, row=row, sticky=tk.W)
            ttk.Entry(self.source_tab, textvariable=var).grid(column=1, row=row)
            row += 1

    def setup_scatterer_tab(self):
        row = 0
        for label, var in self.scatterer_vars.items():
            ttk.Label(self.scatterer_tab, text=label).grid(column=0, row=row, sticky=tk.W)
            ttk.Entry(self.scatterer_tab, textvariable=var).grid(column=1, row=row)
            row += 1

        ttk.Button(self.scatterer_tab, text="Calculate", command=self.update_plot).grid(column=0, row=row, columnspan=2, pady=10)
        ttk.Button(self.scatterer_tab, text="Save as CSV", command=self.save_data_as_csv).grid(column=0, row=row + 1, columnspan=2, pady=10)

    def setup_plot_frame(self):
        self.inputs_frame = ttk.Frame(self.master)
        self.inputs_frame.pack(side=tk.LEFT, fill=tk.Y)
        self.plot_frame = tk.Frame(self.master)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

    def parse_input(self, input_str: str, factor: float = 1) -> float | np.ndarray:
        parts = input_str.split(':')
        if len(parts) == 1:
            return np.array([float(parts[0])]) * factor
        elif len(parts) == 3:
            start, end, points = map(float, parts)
            return np.linspace(start, end, int(points)) * factor
        else:
            raise ValueError("Diameter input must be either a single float or in the format 'start:end:points'")

    def get_data_from_PyMieSim(self) -> None:
        """
        Compute the B1 scattering data using either a single diameter or a range of diameters.
        """
        wavelength_input = self.source_vars['Wavelength (nm)'].get()
        polarization_value_input = self.source_vars['Polarization Value (radians)'].get()
        optical_power_input = self.source_vars['Optical Power (mW)'].get()
        NA_input = self.source_vars['Numerical Aperture (NA)'].get()

        diameter_input = self.scatterer_vars['Diameter (nm)'].get()
        index_input = self.scatterer_vars['Index'].get()
        n_medium_input = self.scatterer_vars['Medium Refractive Index'].get()

        # Process inputs
        wavelength = self.parse_input(input_str=wavelength_input, factor=1e-9)
        polarization_value = self.parse_input(input_str=polarization_value_input, factor=1)
        optical_power = self.parse_input(input_str=optical_power_input, factor=1e-3)
        NA = self.parse_input(input_str=NA_input, factor=1)

        diameters = self.parse_input(input_str=diameter_input, factor=1e-9)
        index = self.parse_input(input_str=index_input, factor=1)
        n_medium = self.parse_input(input_str=n_medium_input, factor=1)

        # Setup PyMieSim components
        self.source_set = Gaussian(
            wavelength=wavelength,
            polarization_value=polarization_value,
            polarization_type='linear',
            optical_power=optical_power,
            NA=NA
        )

        self.scatterer_set = Cylinder(
            diameter=diameters,
            index=index,
            n_medium=n_medium,
            source_set=self.source_set
        )

        self.experiment = Setup(scatterer_set=self.scatterer_set, source_set=self.source_set)

    def save_data_as_csv(self):
        """
        Triggered by the "Save as CSV" button. Opens a file dialog to save the computed data as a CSV file.
        """
        if hasattr(self, 'data'):
            filepath = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
            if filepath:
                # Assuming self.data is a pandas DataFrame or can be converted to one
                self.data.to_csv(filepath, index=False)
                print(f"Data saved to {filepath}")
        else:
            print("No data to save. Please calculate first.")

    def update_plot(self):
        try:
            # Parse inputs from the GUI
            self.get_data_from_PyMieSim()

            # Assuming measure.Qsca can be plotted directly
            self.data = self.experiment.get(measure.Qsca)

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

        except ValueError as e:
            messagebox.showerror("Input Error", str(e))


def main():
    root = tk.Tk()
    PyMieSimGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
