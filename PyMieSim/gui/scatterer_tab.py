
from typing import NoReturn
import tkinter
from PyMieSim.experiment import scatterer
from PyMieSim.gui.base_tab import BaseTab, Widget


class ScattererTab(BaseTab):
    """
    The Scatterer tab for configuring the scatterer in the simulation.

    Inherits from Tab and adds specific controls for configuring the scatterer, such as
    type, diameter, refractive index, and medium refractive index.
    """

    def __init__(self, *args, source_tab, **kwargs):
        """
        Initializes the Scatterer tab with predefined variables for the simulation's scatterer configuration.
        """
        self.type_button = tkinter.StringVar(value='Sphere')

        # diameter_widget = Widget(
        #     default_value='1000',
        #     label='Diameter [nm]',
        #     component_label='diameter',
        #     multiplicative_factor=1e-9
        # )

        # index_widget = Widget(
        #     default_value='1.4',
        #     label='Refractive Index',
        #     component_label='index',
        #     multiplicative_factor=None
        # )

        # n_medium_widget = Widget(
        #     default_value='1.0',
        #     label='Medium Refractive Index',
        #     component_label='n_medium',
        #     multiplicative_factor=None
        # )

        # self.sphere_variables = {
        #     'Diameter [nm]': diameter_widget,
        #     'Index': index_widget,
        #     'Medium Refractive Index': n_medium_widget,
        # }

        self.sphere_variables = {
            'Diameter [nm]': dict(user_input=tkinter.StringVar(value='1000'), factor=1e-9, component_name='diameter'),
            'Index': dict(user_input=tkinter.StringVar(value='1.4'), factor=None, component_name='index'),
            'Medium Refractive Index': dict(user_input=tkinter.StringVar(value='1'), factor=None, component_name='n_medium'),
        }

        self.cylinder_variables = {
            'Diameter [nm]': dict(user_input=tkinter.StringVar(value='1000'), factor=1e-9, component_name='diameter'),
            'Index': dict(user_input=tkinter.StringVar(value='1.4'), factor=None, component_name='index'),
            'Medium Refractive Index': dict(user_input=tkinter.StringVar(value='1'), factor=None, component_name='n_medium'),
        }

        self.coreshell_variables = {
            'Core Diameter [nm]': dict(user_input=tkinter.StringVar(value='1000'), factor=1e-9, component_name='core_diameter'),
            'Shell Width [nm]': dict(user_input=tkinter.StringVar(value='1000'), factor=1e-9, component_name='shell_width'),
            'Core Index': dict(user_input=tkinter.StringVar(value='1.4'), factor=None, component_name='core_index'),
            'Shell Index': dict(user_input=tkinter.StringVar(value='1.5'), factor=None, component_name='shell_index'),
            'Medium Refractive Index': dict(user_input=tkinter.StringVar(value='1'), factor=None, component_name='n_medium'),
        }

        self.source_tab = source_tab

        super().__init__(*args, **kwargs)

        combobox = tkinter.ttk.Combobox(
            self.frame,
            textvariable=self.type_button,
            values=['Sphere', 'Cylinder', 'CoreShell'],
            state="readonly"
        )

        combobox.pack(side=tkinter.TOP)

        combobox.bind("<<ComboboxSelected>>", self.on_type_change)

    def on_type_change(self, event=None) -> NoReturn:
        self.clear_button()

        self.setup_tab()

        self.main_window.axis_tab.clear_button()

        self.main_window.axis_tab.setup_tab()

    @property
    def selected_type(self) -> str:
        return self.type_button.get()

    @property
    def mapping(self) -> str:
        match self.selected_type:
            case 'Sphere':
                return dict(
                    diameter=self.component.diameter,
                    index=self.component.index,
                    n_medium=self.component.n_medium,
                )

            case 'Cylinder':
                return dict(
                    diameter=self.component.diameter,
                    index=self.component.index,
                    n_medium=self.component.n_medium,
                )

            case 'CoreShell':
                return dict(
                    core_diameter=self.component.core_diameter,
                    shell_width=self.component.shell_width,
                    core_index=self.component.core_index,
                    shell_index=self.component.shell_index,
                    n_medium=self.component.n_medium,
                )
            case _:
                raise ValueError('Scatterer type not valid')

    @property
    def variables(self) -> dict:
        match self.selected_type:
            case 'Sphere':
                return self.sphere_variables.copy()
            case 'Cylinder':
                return self.sphere_variables.copy()
            case 'CoreShell':
                return self.coreshell_variables.copy()
            case _:
                raise ValueError('Scatterer type not valid')

    def clear_button(self) -> NoReturn:
        for element in self.non_permanent_widget:
            element.destroy()

    def setup_tab(self):
        """
        Sets up the GUI elements for the Scatterer tab, including labels, entry fields,
        and comboboxes for each scatterer parameter. This method customizes the user interface
        for scatterer configuration, facilitating input of parameters like type, diameter,
        and refractive indices.
        """
        self.non_permanent_widget = []
        for label, var in self.variables.items():
            label = tkinter.Label(self.frame, text=label)
            label.pack(side=tkinter.BOTTOM)
            self.non_permanent_widget.append(label)

            button = tkinter.Entry(self.frame, textvariable=var['user_input'])
            button.pack(side=tkinter.BOTTOM)
            self.non_permanent_widget.append(button)

        self.update_user_input()

        self.setup_component()

    def setup_sphere_component(self):
        self.component = scatterer.Sphere(
            diameter=self.variables['Diameter [nm]']['values'],
            index=self.variables['Index']['values'],
            n_medium=self.variables['Medium Refractive Index']['values'],
            source_set=self.source_tab.component
        )

    def setup_cylinder_component(self):
        self.component = scatterer.Cylinder(
            diameter=self.variables['Diameter [nm]']['values'],
            index=self.variables['Index']['values'],
            n_medium=self.variables['Medium Refractive Index']['values'],
            source_set=self.source_tab.component
        )

    def setup_coreshell_component(self):
        self.component = scatterer.CoreShell(
            core_diameter=self.variables['Core Diameter [nm]']['values'],
            shell_width=self.variables['Shell Width [nm]']['values'],
            core_index=self.variables['Core Index']['values'],
            shell_index=self.variables['Shell Index']['values'],
            n_medium=self.variables['Medium Refractive Index']['values'],
            source_set=self.source_tab.component
        )

    def setup_component(self) -> NoReturn:
        self.update_user_input()

        match self.selected_type:
            case 'Sphere':
                self.setup_sphere_component()
            case 'Cylinder':
                self.setup_cylinder_component()
            case 'CoreShell':
                self.setup_coreshell_component()
            case _:
                raise ValueError('Scatterer type not valid')

