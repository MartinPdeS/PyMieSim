
from typing import NoReturn
import tkinter
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.gui.base_tab import BaseTab


class DetectorTab(BaseTab):
    """
    The Detector tab for configuring the detector in the simulation.
    Inherits from Tab and adds specific controls for configuring the detector, such as
    numerical aperture, gamma, phi, filter type, and coherence.
    """

    def __init__(self, *args, **kwargs):
        self.variables = {
            'NA': dict(user_input=tkinter.StringVar(value='0.2, 0.4'), factor=None),  # Can be ranged
            'Gamma': dict(user_input=tkinter.StringVar(value='0'), factor=None),  # Can be ranged
            'Phi': dict(user_input=tkinter.StringVar(value='0'), factor=None),  # Can be ranged
            'Filter': dict(user_input=tkinter.StringVar(value='0'), factor=None),  # Can be ranged
            'Coherent': dict(user_input=tkinter.BooleanVar(value=True), factor=None, to_float=False),  # True or False
        }

        super().__init__(*args, **kwargs)

    def setup_tab(self):
        """
        Sets up the GUI elements for the Detector tab, including labels, entry fields,
        and checkboxes for each detector parameter. This method enables the configuration
        of detector characteristics such as numerical aperture and coherence.
        """
        for label, var in self.variables.items():
            label = tkinter.Label(self.frame, text=label)
            label.pack(side=tkinter.BOTTOM)

            if isinstance(var['user_input'], tkinter.BooleanVar):  # Special case for Boolean input
                check_button = tkinter.Checkbutton(self.frame, variable=var)
                check_button.pack(side=tkinter.BOTTOM)
            else:
                entry = tkinter.Entry(self.frame, textvariable=var['user_input'])
                entry.pack(side=tkinter.BOTTOM)

        self.setup_component()

    def setup_component(self) -> NoReturn:
        self.update_user_input()

        self.component = Photodiode(
            NA=self.variables['NA']['values'],
            gamma_offset=self.variables['Gamma']['values'],
            phi_offset=self.variables['Phi']['values'],
            polarization_filter=self.variables['Filter']['values'],
            sampling=300
        )

        self.mapping = dict(
            NA=self.component.NA,
            gamma=self.component.gamma_offset,
            phi=self.component.phi_offset,
            polarization_filter=self.component.polarization_filter
        )
