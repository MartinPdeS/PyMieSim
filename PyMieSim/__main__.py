
from typing import NoReturn
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# Assuming PyMieSim and its dependencies are installed
from PyMieSim.experiment.scatterer import Cylinder
from PyMieSim.experiment.source import Gaussian
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.experiment import Setup
from PyMieSim import measure


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
        self.master.title("PyMieSim Graphic Interface")
        self.customize_notebook_style()
        self.setup_variables()
        self.setup_plot_frame()
        self.setup_notebook()
        self.setup_controls()  # Setup controls outside of the tabs

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
            background=[("selected", "#e0e0e0")],
            expand=[("selected", [1, 1, 1, 0])]
        )

        # Style for larger buttons
        style.configure(
            "Large.TButton", 
            font=('Helvetica', 18),    # Font size for buttons
            padding=[20, 15]           # Increased padding for top and bottom to make the button taller
        )

    def setup_variables(self) -> NoReturn:
        """
        Sets up the Tkinter variables for source, scatterer, and detector configurations.
        """
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

        self.detector_vars = {
            'NA': tk.StringVar(value='0.2'),  # Can be ranged
            'Gamma': tk.StringVar(value='0'),  # Can be ranged
            'Phi': tk.StringVar(value='0'),  # Can be ranged
            'Filter': tk.StringVar(value='0'),  # Can be ranged
            'Coherent': tk.BooleanVar(value=True),  # True or False
        }

        self.axis_options = ["Wavelength", "Diameter", "Index", "N_medium", "Gamma", "Phi"]

        self.axis_config_vars = {
            'x axis': tk.StringVar(value=self.axis_options[0]),
            'y axis': tk.StringVar(value=self.axis_options[1]),
        }

    def setup_notebook(self) -> NoReturn:
        """
        Sets up the notebook widget with tabs for Source, Scatterer, and Detector configurations.
        """
        self.notebook = ttk.Notebook(self.master)
        self.notebook.pack(fill=tk.BOTH, expand=True, side=tk.TOP)

        self.source_tab = ttk.Frame(self.notebook)
        self.scatterer_tab = ttk.Frame(self.notebook)
        self.detector_tab = ttk.Frame(self.notebook)
        self.axis_config_tab = ttk.Frame(self.notebook)

        self.notebook.add(self.source_tab, text='Source')
        self.notebook.add(self.scatterer_tab, text='Scatterer')
        self.notebook.add(self.detector_tab, text='Detector')
        self.notebook.add(self.axis_config_tab, text='Axis Configuration')

        self.setup_source_tab()
        self.setup_scatterer_tab()
        self.setup_detector_tab()
        self.setup_axis_config_tab()

    def setup_controls(self) -> NoReturn:
        """
        Sets up control buttons for calculating results and saving data.
        """
        self.controls_frame = ttk.Frame(self.master)
        self.controls_frame.pack(fill=tk.X, side=tk.BOTTOM)

        ttk.Button(
            self.controls_frame,
            text="Calculate",
            style="Large.TButton",  # Apply the custom style for larger buttons
            command=self.update_plot
        ).pack(side=tk.LEFT, padx=5, pady=5)

        ttk.Button(
            self.controls_frame,
            text="Save as CSV",
            style="Large.TButton",  # Apply the custom style for larger buttons
            command=self.save_data_as_csv
        ).pack(side=tk.LEFT, padx=5, pady=5)

    def setup_axis_config_tab(self):
        ttk.Label(
            self.axis_config_tab,
            text="x axis"
        ).grid(column=0, row=0, sticky=tk.W, padx=5, pady=5)

        ttk.Combobox(
            self.axis_config_tab,
            textvariable=self.axis_config_vars['x axis'],
            values=self.axis_options,
            state="readonly"
        ).grid(column=1, row=0, padx=5, pady=5)

        ttk.Label(
            self.axis_config_tab,
            text="y axis"
        ).grid(column=0, row=1, sticky=tk.W, padx=5, pady=5)

        ttk.Combobox(
            self.axis_config_tab,
            textvariable=self.axis_config_vars['y axis'],
            values=self.axis_options,
            state="readonly"
        ).grid(column=1, row=1, padx=5, pady=5)

    def setup_source_tab(self) -> NoReturn:
        """
        Sets up widgets for the Source configuration tab.
        """
        row = 0
        for label, var in self.source_vars.items():
            ttk.Label(self.source_tab, text=label).grid(column=0, row=row, sticky=tk.W)
            ttk.Entry(self.source_tab, textvariable=var).grid(column=1, row=row)
            row += 1

    def setup_scatterer_tab(self) -> NoReturn:
        """
        Sets up widgets for the Scatterer configuration tab.
        """
        row = 0
        for label, var in self.scatterer_vars.items():
            ttk.Label(self.scatterer_tab, text=label).grid(column=0, row=row, sticky=tk.W)
            ttk.Entry(self.scatterer_tab, textvariable=var).grid(column=1, row=row)
            row += 1

    def setup_detector_tab(self) -> NoReturn:
        """
        Sets up widgets for the Detector configuration tab, including a special case for the Coherent boolean input.
        """
        row = 0
        for label, var in self.detector_vars.items():
            ttk.Label(self.detector_tab, text=label).grid(column=0, row=row, sticky=tk.W)

            if label == 'Coherent':  # Special case for Boolean input
                ttk.Checkbutton(
                    self.detector_tab,
                    text='',
                    variable=var,
                    onvalue=True,
                    offvalue=False
                ).grid(column=1, row=row)
            else:
                ttk.Entry(self.detector_tab, textvariable=var).grid(column=1, row=row)
            row += 1

    def setup_plot_frame(self) -> NoReturn:
        """
        Sets up the frame for displaying plots.
        """
        self.plot_frame = tk.Frame(self.master)
        self.plot_frame.pack(fill=tk.BOTH, expand=True)

    def parse_input(self, input_str: str, factor: float = 1) -> float | np.ndarray:
        """
        Parses input strings to handle single values or ranges. Supports notation for ranges (e.g., "start:end:points").

        Parameters:
            input_str (str): The input string to parse.
            factor (float): Multiplication factor for unit conversion.

        Returns:
            float | np.ndarray: The parsed input as a single float or a numpy array.
        """
        if ":" in input_str:
            parts = input_str.split(':')
            start, end, points = map(float, parts)
            return np.linspace(start, end, int(points)) * factor

        if "," in input_str:
            parts = input_str.split(',')
            values = [float(value) for value in parts]
            return np.asarray(values) * factor

        return np.asarray(float(input_str)) * factor

    def set_source(self) -> NoReturn:
        wavelength_input = self.source_vars['Wavelength (nm)'].get()
        polarization_value_input = self.source_vars['Polarization Value (radians)'].get()
        optical_power_input = self.source_vars['Optical Power (mW)'].get()
        NA_input = self.source_vars['Numerical Aperture (NA)'].get()

        # Process inputs
        wavelength = self.parse_input(input_str=wavelength_input, factor=1e-9)
        polarization_value = self.parse_input(input_str=polarization_value_input, factor=1)
        optical_power = self.parse_input(input_str=optical_power_input, factor=1e-3)
        NA = self.parse_input(input_str=NA_input, factor=1)

        self.source_set = Gaussian(
            wavelength=wavelength,
            polarization_value=polarization_value,
            polarization_type='linear',
            optical_power=optical_power,
            NA=NA
        )

    def set_scatterer(self) -> NoReturn:
        diameter_input = self.scatterer_vars['Diameter (nm)'].get()
        index_input = self.scatterer_vars['Index'].get()
        n_medium_input = self.scatterer_vars['Medium Refractive Index'].get()

        diameters = self.parse_input(input_str=diameter_input, factor=1e-9)
        index = self.parse_input(input_str=index_input, factor=1)
        n_medium = self.parse_input(input_str=n_medium_input, factor=1)

        self.scatterer_set = Cylinder(
            diameter=diameters,
            index=index,
            n_medium=n_medium,
            source_set=self.source_set
        )

    def set_detector(self) -> NoReturn:
        NA_input = self.detector_vars['NA'].get()
        gamma_input = self.detector_vars['Gamma'].get()
        phi_input = self.detector_vars['Phi'].get()
        filter_input = self.detector_vars['Filter'].get()

        NA = self.parse_input(input_str=NA_input, factor=1)
        gamma = self.parse_input(input_str=gamma_input, factor=1)
        phi = self.parse_input(input_str=phi_input, factor=1)
        polarization_filter = self.parse_input(input_str=filter_input, factor=1)

        self.detector_set = Photodiode(
            NA=NA,
            gamma_offset=gamma,
            phi_offset=phi,
            polarization_filter=polarization_filter,
            sampling=300
        )

    def get_data_from_PyMieSim(self) -> NoReturn:
        """
        Compute the B1 scattering data using either a single diameter or a range of diameters.
        """
        self.set_source()

        self.set_scatterer()

        self.set_detector()

        self.axis_mapping = {
            'Diameter': self.scatterer_set.diameter,
            'Index': self.scatterer_set.index,
            'N_medium': self.scatterer_set.n_medium,
            'Wavelength': self.source_set.wavelength,
            "Gamma": self.detector_set.gamma_offset,
            "Phi": self.detector_set.phi_offset
        }

        self.experiment = Setup(
            scatterer_set=self.scatterer_set,
            source_set=self.source_set,
            detector_set=self.detector_set
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

    def update_plot(self) -> NoReturn:
        try:
            # Parse inputs from the GUI
            self.get_data_from_PyMieSim()

            # Assuming measure.Qsca can be plotted directly
            self.data = self.experiment.get(measure.coupling)

            x_axis = self.axis_config_vars['x axis'].get()
            x_axis = self.axis_mapping[x_axis]

            figure = self.data.plot(x=x_axis)
            figure.unit_size = (7, 4)
            figure._render_()
            figure = figure._mpl_figure

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
