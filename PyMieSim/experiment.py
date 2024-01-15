#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from dataclasses import dataclass

from DataVisual import Xparameter
from DataVisual import DataV
import DataVisual.tables as Table

from PyMieSim import load_lp_mode
from PyMieSim.physics import power_to_amplitude
from PyMieSim.polarization import LinearPolarization
from PyMieSim.binary.Experiment import CppExperiment
import PyMieSim.datavisual_x_parameters as Kwargs

from PyMieSim.binary.Sets import (
    CppCoreShellSet,
    CppCylinderSet,
    CppSphereSet,
    CppSourceSet,
    CppDetectorSet
)


@dataclass
class CoreShellSet():
    core_diameter: list
    """ diameter of the core of the single scatterer [m]. """
    shell_width: list
    """ Width of the shell of the single scatterer [m]. """
    core_index: tuple = tuple()
    """ Refractive index of the core of the scatterer. """
    shell_index: tuple = tuple()
    """ Refractive index of the shell of the scatterer. """
    core_material: tuple = tuple()
    """ Core material of which the scatterer is made of. Only if core_index is not specified.  """
    shell_material: tuple = tuple()
    """ Shell material of which the scatterer is made of. Only if shell_index is not specified.  """
    n_medium: list = 1.0
    """ Refractive index of scatterer medium. """
    name: str = 'coreshell'
    """name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        super().__init__()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.n_medium = Xparameter(values=self.n_medium, **Kwargs.n_medium)

        self.bounded_core = True if len(self.core_material) != 0 else False

        self.bounded_shell = True if len(self.shell_material) != 0 else False

        self.core_diameter = Xparameter(values=self.core_diameter, **Kwargs.core_diameter)

        self.shell_width = Xparameter(values=self.shell_width, **Kwargs.shell_width)

        self.core_material = Xparameter(values=self.core_material, **Kwargs.core_material)

        self.shell_material = Xparameter(values=self.shell_material, **Kwargs.shell_material)

        self.core_index = Xparameter(values=self.core_index, **Kwargs.core_index)

        self.shell_index = Xparameter(values=self.shell_index, **Kwargs.shell_index)

    def bind_to_experiment(self, experiment):
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_coreshell(self.binding)

    def get_core_shell_diameter_from_shell_width(
            self,
            core_diameter: numpy.ndarray,
            shell_width: numpy.ndarray) -> tuple:
        """
        Gets the core shell diameter from shell width.

        :param      core_diameter:  The core diameter
        :type       core_diameter:  numpy.ndarray
        :param      shell_width:    The shell width
        :type       shell_width:    numpy.ndarray

        :returns:   The core shell diameter from shell width.
        :rtype:     tuple
        """
        _core_diameter = numpy.tile(core_diameter, reps=len(shell_width))
        shell_diameter = numpy.add.outer(core_diameter, shell_width)
        return _core_diameter.flatten(), shell_diameter.T.flatten()

    def bind_core_material_shell_material(self, source) -> None:
        core_material = [
            material.GetRI(source.values) for material in self.core_material
        ]

        shell_material = [
            material.GetRI(source.values) for material in self.shell_material
        ]

        core_material = numpy.asarray(core_material).astype(complex)
        shell_material = numpy.asarray(shell_material).astype(complex)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_material=core_material,
            shell_material=shell_material,
            n_medium=self.n_medium.values,
        )

    def bind_core_material_shell_index(self, source) -> None:
        core_material = [
            material.GetRI(source.values) for material in self.core_material
        ]

        core_material = numpy.asarray(core_material).astype(complex)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_material=core_material,
            shell_index=self.shell_index.values,
            n_medium=self.n_medium.values,
        )

    def bind_index_material_shell_material(self, source) -> None:
        shell_material = [
            material.GetRI(source.values) for material in self.shell_material
        ]

        shell_material = numpy.asarray(shell_material).astype(complex)

        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_index=self.core_index.values,
            shell_material=shell_material,
            n_medium=self.n_medium.values,
        )

    def bind_index_material_shell_index(self, source) -> None:
        self.binding = CppCoreShellSet(
            core_diameter=self.core_diameter.values,
            shell_width=self.shell_width.values,
            core_index=self.core_index.values,
            shell_index=self.shell_index.values,
            n_medium=self.n_medium.values,
        )

    def evaluate_index_material(self, source) -> None:
        if self.bounded_core and self.bounded_shell:
            self.bind_core_material_shell_material(source)

        if self.bounded_core and not self.bounded_shell:
            self.bind_core_material_shell_index(source)

        if not self.bounded_core and self.bounded_shell:
            self.bind_index_material_shell_material(source)

        if not self.bounded_core and not self.bounded_shell:
            self.bind_index_material_shell_index(source)

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.core_diameter = numpy.atleast_1d(self.core_diameter).astype(float)
        self.shell_width = numpy.atleast_1d(self.shell_width).astype(float)

        self.core_index = numpy.atleast_1d(self.core_index).astype(numpy.complex128)
        self.shell_index = numpy.atleast_1d(self.shell_index).astype(numpy.complex128)

        self.n_medium = numpy.atleast_1d(self.n_medium).astype(float)

        self.core_material = numpy.atleast_1d(self.core_material)
        self.shell_material = numpy.atleast_1d(self.shell_material)

    def append_to_table(self, table) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        if self.bounded_core and self.bounded_shell:
            index_material_table = [self.core_material, self.shell_material]

        if self.bounded_core and not self.bounded_shell:
            index_material_table = [self.core_material, self.shell_index]

        if not self.bounded_core and self.bounded_shell:
            index_material_table = [self.core_index, self.shell_material]

        if not self.bounded_core and not self.bounded_shell:
            index_material_table = [self.core_index, self.shell_index]

        return [*table, self.core_diameter, self.shell_width, *index_material_table, self.n_medium]


@dataclass
class SphereSet():
    diameter: list
    """ diameter of the single scatterer in unit of meter. """
    index: list = tuple()
    """ Refractive index of scatterer. """
    material: list = tuple()
    """ material of which the scatterer is made of. Only if index is not specified. """
    n_medium: list = 1.0
    """ Refractive index of scatterer medium. """
    name: str = 'sphere'
    """name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        super().__init__()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.n_medium = Xparameter(values=self.n_medium, **Kwargs.n_medium)

        self.bounded_index = True if len(self.material) != 0 else False

        self.diameter = Xparameter(values=self.diameter, **Kwargs.diameter)

        self.material = Xparameter(values=self.material, **Kwargs.material)

        self.index = Xparameter(values=self.index, **Kwargs.index)

    def bind_to_experiment(self, experiment) -> None:
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_sphere(self.binding)

    def evaluate_index_material(self, Source) -> None:
        if self.bounded_index:
            material = [
                material.GetRI(Source.values) for material in self.material
            ]

            material_index = numpy.asarray(material).astype(complex)

            self.binding = CppSphereSet(
                diameter=self.diameter.values.astype(float),
                material=material_index,
                n_medium=self.n_medium.values.astype(float)
            )

        else:
            self.binding = CppSphereSet(
                diameter=self.diameter.values,
                index=self.index.values,
                n_medium=self.n_medium.values
            )

    def append_to_table(self, table) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        if self.bounded_index:
            index_material_table = [self.material]

        else:
            index_material_table = [self.index]

        return [*table, self.diameter, *index_material_table, self.n_medium]

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.diameter = numpy.atleast_1d(self.diameter).astype(float)

        self.index = numpy.atleast_1d(self.index).astype(numpy.complex128)

        self.n_medium = numpy.atleast_1d(self.n_medium).astype(float)

        self.material = numpy.atleast_1d(self.material)


@dataclass
class CylinderSet():
    diameter: list
    """ diameter of the single scatterer in unit of meter. """
    index: tuple = tuple()
    """ Refractive index of scatterer. """
    material: tuple = tuple()
    """ Refractive index of scatterer medium. """
    n_medium: list = 1.0
    """ material of which the scatterer is made of. Only if index is not specified. """
    name: str = 'cylinder'
    """name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        super().__init__()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.n_medium = Xparameter(values=self.n_medium, **Kwargs.n_medium)

        self.bounded_index = True if len(self.material) != 0 else False

        self.diameter = Xparameter(values=self.diameter, **Kwargs.diameter)

        self.material = Xparameter(values=self.material, **Kwargs.material)

        self.index = Xparameter(values=self.index, **Kwargs.index)

    def bind_to_experiment(self, experiment):
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_cylinder(self.binding)

    def evaluate_index_material(self, source):
        if self.bounded_index:
            material = numpy.asarray([material.GetRI(source.values) for material in self.material])

            self.binding = CppCylinderSet(
                diameter=self.diameter.values,
                material=material,
                n_medium=self.n_medium.values
            )

        else:
            self.binding = CppCylinderSet(
                diameter=self.diameter.values,
                index=self.index.values,
                n_medium=self.n_medium.values
            )

    def append_to_table(self, table):
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        if self.bounded_index:
            index_material_table = [self.material]

        else:
            index_material_table = [self.index]

        return [*table, self.diameter, *index_material_table, self.n_medium]

    def format_inputs(self):
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.diameter = numpy.atleast_1d(self.diameter).astype(float)

        self.index = numpy.atleast_1d(self.index).astype(numpy.complex128)

        self.n_medium = numpy.atleast_1d(self.n_medium).astype(float)

        self.material = numpy.atleast_1d(self.material)


@dataclass
class SourceSet(object):
    wavelength: float
    """ Wavelenght of the light field. """
    optical_power: float
    """ Optical power in unit of Watt """
    NA: float
    """ Numerical aperture of the source """
    linear_polarization: float = None
    """ polarization of the light field in degree. """
    amplitude: float = None
    """ Maximal value of the electric field at focus point. """
    name: str = 'PlaneWave'
    """ name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.amplitude = power_to_amplitude(
            wavelength=self.wavelength,
            optical_power=self.optical_power,
            NA=self.NA,
        )

        self.build_x_parameters()

        self.generate_binding()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.wavelength = Xparameter(values=self.wavelength, **Kwargs.wavelength)

        self.linear_polarization = Xparameter(
            values=self.jones_vector,
            representation=self.linear_polarization.angle_list,
            **Kwargs.linear_polarization
        )

    def generate_binding(self) -> None:
        """
        Generate the C++ binding

        :returns:   No return
        :rtype:     None
        """
        self.binding = CppSourceSet(
            wavelength=self.wavelength.values,
            jones_vector=self.jones_vector,
            amplitude=self.amplitude
        )

    def format_inputs(self):
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        if numpy.iterable(self.linear_polarization):
            self.linear_polarization = LinearPolarization(*self.linear_polarization)
        else:
            self.linear_polarization = LinearPolarization(self.linear_polarization)

        self.wavelength = numpy.atleast_1d(self.wavelength).astype(float)

        self.amplitude = numpy.atleast_1d(self.amplitude).astype(float)

        self.jones_vector = numpy.atleast_1d(self.linear_polarization.jones_vector).astype(numpy.complex128).T

    def bind_to_experiment(self, experiment):
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_source(self.binding)

    def append_to_table(self, table):
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, self.wavelength, self.linear_polarization]


@dataclass
class PhotodiodeSet():
    NA: list
    """ Numerical aperture of imaging system. """
    gamma_offset: list
    """ Angle [degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: list
    """ Angle [degree] offset of detector in the direction parallel to polarization. """
    polarization_filter: list
    """ Angle [degree] of polarization filter in front of detector. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = False
    """ Describe the detection scheme coherent or uncoherent. """
    sampling: int = 200
    """ Sampling number for the field evaluation. """
    name = "Photodiode"
    """ name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        self.initialize_binding()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        self.scalarfield = numpy.ones([1, self.sampling])

        self.NA = Xparameter(values=self.NA, **Kwargs.NA)

        self.phi_offset = Xparameter(values=self.phi_offset, **Kwargs.phi_offset)

        self.gamma_offset = Xparameter(values=self.gamma_offset, **Kwargs.gamma_offset)

        self.polarization_filter = Xparameter(values=self.polarization_filter, **Kwargs.polarization_filter)

        self.scalarfield = Xparameter(
            values=self.scalarfield,
            representation=numpy.asarray(['Photodiode']),
            **Kwargs.scalarfield
        )

    def initialize_binding(self):
        point_coupling = True if self.coupling_mode == 'point' else False
        phi_offset_rad = numpy.deg2rad(self.phi_offset.values)
        gamma_offset_rad = numpy.deg2rad(self.gamma_offset.values)

        self.binding = CppDetectorSet(
            scalarfield=self.scalarfield.values,
            NA=self.NA.values,
            phi_offset=phi_offset_rad,
            gamma_offset=gamma_offset_rad,
            polarization_filter=self.polarization_filter.values,
            point_coupling=point_coupling,
            coherent=self.coherent
        )

    def bind_to_experiment(self, experiment):
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_detector(self.binding)

    def append_to_table(self, table):
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, self.scalarfield, self.NA, self.phi_offset, self.gamma_offset, self.polarization_filter]

    def format_inputs(self):
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.NA = numpy.atleast_1d(self.NA).astype(float)

        self.phi_offset = numpy.atleast_1d(self.phi_offset).astype(float)

        self.gamma_offset = numpy.atleast_1d(self.gamma_offset).astype(float)

        self.polarization_filter = numpy.atleast_1d(self.polarization_filter).astype(float)


@dataclass
class LPModeSet():
    mode_number: list
    """ List of mode to be used. """
    NA: list
    """ Numerical aperture of imaging system. """
    gamma_offset: list
    """ Angle [degree] offset of detector in the direction perpendicular to polarization. """
    phi_offset: list
    """ Angle [degree] offset of detector in the direction parallel to polarization. """
    polarization_filter: list
    """ Angle [degree] of polarization filter in front of detector. """
    coupling_mode: str = 'point'
    """ Method for computing mode coupling. Either point or Mean. """
    coherent: bool = True
    """ Describe the detection scheme coherent or uncoherent. """
    sampling: int = 200
    """ Sampling number for the field evaluation. """
    name = "LPMode"
    """ name of the set """

    def __post_init__(self):
        self.format_inputs()

        self.build_x_parameters()

        self.initialize_binding()

    def build_x_parameters(self) -> None:
        """
        Builds the parameters that will be passed in XTable for DataVisual.

        :returns:   No Return
        :rtype:     None
        """
        representation = [
            f"{mode}" for mode in self.mode_number
        ]

        self.NA = Xparameter(values=self.NA, **Kwargs.NA)

        self.phi_offset = Xparameter(values=self.phi_offset, **Kwargs.phi_offset)

        self.gamma_offset = Xparameter(values=self.gamma_offset, **Kwargs.gamma_offset)

        self.polarization_filter = Xparameter(values=self.polarization_filter, **Kwargs.polarization_filter)

        self.scalarfield = Xparameter(
            values=self.get_fields_array(),
            representation=representation,
            **Kwargs.scalarfield
        )

    def get_fields_array(self) -> numpy.ndarray:
        """
        Gets the field arrays for the detection schemes the first dimension is the different fields.
        Second dimension are the individual mode fields.

        :returns:   The fields array.
        :rtype:     numpy.ndarray
        """
        fields_array = [
            load_lp_mode(
                mode_number=mode,
                sampling=self.sampling,
                structure_type='unstructured'
            ) for mode in self.mode_number
        ]

        return numpy.asarray(fields_array).astype(numpy.complex128)

    def initialize_binding(self) -> None:
        """
        Initializes the cpp binding of the LPmode detector set.

        :returns:   No return
        :rtype:     None
        """
        point_coupling = True if self.coupling_mode == 'point' else False
        phi_offset_rad = numpy.deg2rad(self.phi_offset.values)
        gamma_offset_rad = numpy.deg2rad(self.gamma_offset.values)

        self.binding = CppDetectorSet(
            scalarfield=self.scalarfield.values,
            NA=self.NA.values,
            phi_offset=phi_offset_rad,
            gamma_offset=gamma_offset_rad,
            polarization_filter=self.polarization_filter.values,
            point_coupling=point_coupling,
            coherent=self.coherent
        )

    def bind_to_experiment(self, experiment) -> None:
        """
        Bind this specific set to a Setup experiment.

        :param      experiment:  The experiment
        :type       experiment:  Setup

        :returns:   No return
        :rtype:     None
        """
        experiment.binding.set_detector(self.binding)

    def append_to_table(self, table: list) -> list:
        """
        Append elements to the xTable from the DataVisual library for the plottings.

        :param      table:  The table
        :type       table:  list

        :returns:   The updated list
        :rtype:     list
        """
        return [*table, self.scalarfield, self.NA, self.phi_offset, self.gamma_offset, self.polarization_filter]

    def format_inputs(self) -> None:
        """
        Format the inputs given by the user into numpy array. Those inputs are subsequently
        sent to the cpp binding.

        :returns:   No return
        :rtype:     None
        """
        self.mode_number = numpy.atleast_1d(self.mode_number).astype(str)

        self.NA = numpy.atleast_1d(self.NA).astype(float)

        self.phi_offset = numpy.atleast_1d(self.phi_offset).astype(float)

        self.gamma_offset = numpy.atleast_1d(self.gamma_offset).astype(float)

        self.polarization_filter = numpy.atleast_1d(self.polarization_filter).astype(float)


@dataclass
class Setup(object):
    scatterer_set: SphereSet | CylinderSet | CoreShellSet
    """ Scatterer set instance which defined the ranging paremters to be measured """
    source_set: SourceSet
    """ Source set instance which defined the ranging paremters to be measured, source is PlaneWave """
    detector_set: PhotodiodeSet | LPModeSet = None
    """ Detector set instance which defined the ranging paremters to be measured, default is None as all parameter do not need a detector """

    def __post_init__(self):
        self.scatterer_set.evaluate_index_material(self.source_set.wavelength)

        self.binding = CppExperiment()

        self.bind_sets_to_experiment()

        self.x_table = self.source_set.append_to_table(table=[])
        self.x_table = self.scatterer_set.append_to_table(table=self.x_table)

        if self.detector_set is not None:
            self.x_table = self.detector_set.append_to_table(table=self.x_table)

    def bind_sets_to_experiment(self):
        self.source_set.bind_to_experiment(self)
        self.scatterer_set.bind_to_experiment(self)

        if self.detector_set:
            self.detector_set.bind_to_experiment(self)

    def Get(self, measure) -> DataV:
        """
        Compute the measure provided and return a DataV structured array.

        :param      measures:  The measures
        :type       measures:  list

        :returns:   The data structure.
        :rtype:     DataV
        """
        self.y_table = [measure]

        measure_string = f'get_{self.scatterer_set.name}_{measure.name}'

        array = getattr(self.binding, measure_string)()

        x_table = Table.Xtable(self.x_table)

        return DataV(
            array=array,
            x_table=x_table,
            y_parameter=measure
        )


# -
