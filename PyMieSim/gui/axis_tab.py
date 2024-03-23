
from typing import NoReturn
import tkinter
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab


class AxisTab(BaseTab):
    measure_map = scatterer.Sphere.available_measure_list

    def __init__(self, *args, other_tabs: list, **kwargs):
        """
        Initializes the Axis Configuration tab with predefined variables for selecting
        the axes for plotting the simulation results.
        """
        self.other_tabs = other_tabs

        super().__init__(*args, **kwargs)

    def clear_button(self) -> NoReturn:
        for element in self.non_permanent_widget:
            element.destroy()

    def setup_tab(self) -> NoReturn:
        """
        Sets up the GUI elements for the Axis Configuration tab, including labels and comboboxes
        for selecting the x and y axis variables. This method provides a user interface for
        specifying the variables to be plotted on the axes.
        """
        self.x_axis_options = list(self.axis_mapping.keys())
        self.y_axis_options = list(self.measure_map.keys())

        self.variables = {
            'x axis': dict(user_input=tkinter.StringVar(value="wavelength"), factor=None),
            'y axis': dict(user_input=tkinter.StringVar(value='coupling'), factor=None),
        }

        self.non_permanent_widget = []

        label = tkinter.Label(
            self.frame,
            text="x axis"
        )

        label.pack(side=tkinter.BOTTOM)

        self.non_permanent_widget.append(label)

        combox = tkinter.ttk.Combobox(
            self.frame,
            textvariable=self.variables['x axis']['user_input'],
            values=self.x_axis_options,
            state="readonly"
        )

        combox.pack(side=tkinter.BOTTOM)

        self.non_permanent_widget.append(combox)

        label = tkinter.Label(
            self.frame,
            text="y axis"
        )

        label.pack(side=tkinter.BOTTOM)
        self.non_permanent_widget.append(label)

        combox = tkinter.ttk.Combobox(
            self.frame,
            textvariable=self.variables['y axis']['user_input'],
            values=self.y_axis_options,
            state="readonly"
        )

        combox.pack(side=tkinter.BOTTOM)
        self.non_permanent_widget.append(combox)

    @property
    def x_axis(self) -> str:
        x_axis = self.variables['x axis']['user_input'].get()

        return self.axis_mapping[x_axis]

    @property
    def axis_mapping(self):
        _axis_mapping = {}
        for tab in self.other_tabs:
            _axis_mapping.update(tab.mapping)

        return _axis_mapping
