#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.material.base_class import BaseMaterial
from pydantic import field_validator
from pint_pandas import PintArray
import numpy

from TypedUnit import Length, RefractiveIndex
from PyMieSim.binary.interface_sets import CppScattererProperties, CppMediumProperties


class BaseScatterer():
    """
    Base class for scatterer objects.  This class handles the initialization and setup of
    scatterer parameters for use in PyMieSim simulations.

    """
    mapping = None
    binding_kwargs = None
    binding = None

    @field_validator('property', 'core_property', 'shell_property', 'medium_property', mode='plain')
    def _validate_array(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""
        value = numpy.atleast_1d(value)

        if isinstance(value, RefractiveIndex):
            RefractiveIndex.check(value)
            return value

        assert numpy.all([isinstance(p, BaseMaterial) for p in value]), f"{value} must be either a refractive index quantity (RIU) or BaseMaterial."

        return value

    @field_validator('diameter', 'core_diameter', 'shell_thickness', mode='plain')
    def _validate_length_quantity(cls, value):
        """
        Ensures that diameter is Quantity objects with length units.
        """
        Length.check(value)
        return numpy.atleast_1d(value)

    def _add_properties(self, name: str, properties: RefractiveIndex | BaseMaterial) -> None:
        """
        Determines whether the provided property is a refractive index (Quantity) or a material (BaseMaterial),
        and returns the corresponding values.

        Parameters
        ----------
        property : Quantity or BaseMaterial
            The core property to be assigned, which can either be a refractive index (Quantity) or a material (BaseMaterial).

        Raises
        ------
        ValueError:
            If the provided property is neither a Quantity (refractive index) nor a BaseMaterial.
        """
        CPPClass = CppMediumProperties if name == 'medium' else CppScattererProperties

        if all(isinstance(item, RefractiveIndex) for item in properties):
            instance = CPPClass(index_properties=properties.magnitude)
            self.binding_kwargs[f'{name}_properties'] = instance

            return instance

        elif all(isinstance(item, BaseMaterial) for item in properties):
            eval_index = numpy.asarray([m.compute_refractive_index(self.source.wavelength.to('meter').magnitude) for m in properties])
            instance = CPPClass(material_properties=eval_index)
            self.binding_kwargs[f'{name}_properties'] = instance

            return instance

        else:
            raise TypeError("All elements in the list must be of type 'Quantity' or 'BaseMaterial', not a mix of both.")

    def _generate_mapping(self) -> None:
        """
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns
        -------
        list
            A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """
        for attr in self.attributes:
            values = getattr(self, attr)
            if values is None:
                continue

            if hasattr(values, 'magnitude'):
                self.mapping["scatterer:" + attr] = PintArray(values.magnitude, dtype=values.units)
            else:
                self.mapping["scatterer:" + attr] = [repr(m) for m in values]
