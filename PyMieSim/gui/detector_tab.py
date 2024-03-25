#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn
from PyMieSim.experiment.detector import Photodiode
from PyMieSim.gui.base_tab import BaseTab, Widget, WidgetCollection


class DetectorTab(BaseTab):
    """
    The Detector tab for configuring the detector in the simulation.
    Inherits from Tab and adds specific controls for configuring the detector, such as
    numerical aperture, gamma, phi, filter type, and coherence.
    """

    def __init__(self, *args, **kwargs):

        na_widget = Widget(
            default_value='0.2, 0.4',
            label='Numerical aperture (NA)',
            component_label='NA',
        )

        gamma_widget = Widget(
            default_value='0',
            label='Gamma [degree]',
            component_label='gamma_offset',
        )

        phi_widget = Widget(
            default_value='0',
            label='Phi [degree]',
            component_label='phi_offset',
        )

        filter_widget = Widget(
            default_value='0.',
            label='Polarization filter [degree]',
            component_label='polarization_filter',
        )

        self.variables = WidgetCollection(
            na_widget,
            gamma_widget,
            phi_widget,
            filter_widget,
        )

        super().__init__(*args, **kwargs)

    def setup_tab(self):
        """
        Sets up the GUI elements for the Detector tab, including labels, entry fields,
        and checkboxes for each detector parameter. This method enables the configuration
        of detector characteristics such as numerical aperture and coherence.
        """
        self.variables.setup_widgets(frame=self.frame)

        self.setup_component()

    def setup_component(self) -> NoReturn:
        self.update_user_input()

        self.component = Photodiode(**self.variables.to_component_dict(), sampling=300)

        self.mapping = dict(
            NA=self.component.NA,
            gamma=self.component.gamma_offset,
            phi=self.component.phi_offset,
            polarization_filter=self.component.polarization_filter
        )
