
from typing import NoReturn
import tkinter as tk
from tkinter import ttk
from PyMieSim.experiment.source import Gaussian
from PyMieSim.gui.base_tab import BaseTab


class SourceTab(BaseTab):
    """
    The Source tab for configuring the light source in the simulation.

    Inherits from Tab and adds specific controls for configuring the source, such as
    wavelength, polarization, optical power, and numerical aperture.
    """

    def __init__(self, *args, **kwargs):
        """
        Initializes the Source tab with predefined variables for the simulation's source configuration.
        """
        self.variables = {
            'Wavelength [nm]': dict(user_input=tk.StringVar(value='400:950:200'), factor=1e-9),
            'Polarization Value (radians)': dict(user_input=tk.StringVar(value='0'), factor=None),
            'Optical Power (mW)': dict(user_input=tk.StringVar(value='1'), factor=1e-3),
            'Numerical Aperture (NA)': dict(user_input=tk.StringVar(value='0.2'), factor=None),
        }

        super().__init__(*args, **kwargs)

    def setup_tab(self) -> NoReturn:
        """
        Sets up the GUI elements for the Source tab, including labels and entry fields for each source parameter.
        """
        self.update_user_input()

        for label, var in self.variables.items():
            label = ttk.Label(self.frame, text=label)
            label.pack(side=tk.BOTTOM)

            entry = ttk.Entry(self.frame, textvariable=var['user_input'])
            entry.pack(side=tk.BOTTOM)

        self.setup_component()

    def setup_component(self) -> NoReturn:
        self.update_user_input()

        self.component = Gaussian(
            wavelength=self.variables['Wavelength [nm]']['values'],
            polarization_value=self.variables['Polarization Value (radians)']['values'],
            polarization_type='linear',
            optical_power=self.variables['Optical Power (mW)']['values'],
            NA=self.variables['Numerical Aperture (NA)']['values'],
        )

        self.mapping = dict(
            wavelength=self.component.wavelength,
            polarization=self.component.polarization_value
        )
