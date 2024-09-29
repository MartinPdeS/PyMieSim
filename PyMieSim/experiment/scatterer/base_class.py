#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydantic.dataclasses import dataclass
from dataclasses import fields
from pydantic import ConfigDict
import pint_pandas

import numpy
from PyMieSim.experiment import parameters

config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


@dataclass(config=config_dict)
class BaseScatterer():
    """
    Base class for scatterer objects. This class handles the initialization and setup of
    scatterer parameters for use in PyMieSim simulations.

    """
    mapping = None
    binding_kwargs = None
    binding = None

    def _add_material_index_to_mapping(self, name: str, indexes: numpy.ndarray, materials: numpy.ndarray, data_type: type = object) -> None:
        """
        Adds material or refractive index details to a mapping dictionary.

        This method is used to create a mapping of material properties to human-readable and accessible formats
        for UI or data outputs. The key in the mapping dictionary is constructed using the provided name.

        Parameters:
            name (str): The base name to use for the keys in the mapping dictionary. This name is used to differentiate between different materials or indices.
            indexes (numpy.ndarray): The array of refractive index values.
            materials (numpy.ndarray): The array of material objects.
            data_type (type): The expected data type of the material or index values. Default is `object`.

        """
        if materials is not None:
            key = f"{name}_material" if name else "material"
            base_values = materials
        else:
            key = f"{name}_index" if name else "index"
            base_values = indexes

        res = getattr(parameters, key)
        res.base_values = res.values = base_values
        self.mapping[key] = res

    def add_material_index_to_binding_kwargs(self, name: str, indexes: numpy.ndarray, materials: numpy.ndarray, data_type: type = object) -> None:
        """
        Adds either material properties or a refractive index to the binding keyword arguments for the experiment.

        This method validates and processes the material or index information, converting it into a format suitable
        for simulation use, and ensuring that either a material or an index is provided but not both.

        Parameters:
            name (str): The base name for the material or index. This name helps identify the property and is used to handle multiple materials or indices.
            indexes (numpy.ndarray): The array of refractive index values.
            materials (numpy.ndarray): The array of material objects.
            data_type (type): The expected data type of the material or index values. Default is `object`.

        Raises:
            ValueError: If both a material and an index are provided, or if neither is provided.

        Returns:
            None
        """
        if (materials is not None) == (indexes is not None):
            raise ValueError(f"Either {name} material or {name} index must be provided.")

        if materials:
            key = f"{name}_material" if name else "material"
            wavelength = self.source.wavelength
            indexes = numpy.asarray(
                [m.compute_refractive_index(wavelength.to_base_units().magnitude) for m in materials]
            )
        else:
            key = f"{name}_index" if name else "index"
            indexes = indexes

        indexes = numpy.atleast_1d(indexes)

        if data_type == float and numpy.any(numpy.iscomplex(indexes)):
            indexes = indexes.real

        self.binding_kwargs[key] = indexes.astype(data_type)

    def _generate_mapping(self) -> None:
        """
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns:
            list: A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """

        for attr in [f.name for f in fields(self) if f.name != 'source']:
            values = getattr(self, attr)

            if values is None: continue

            # attr = attr.replace('_', ' ').capitalize()

            if hasattr(values, 'magnitude'):
                magnitude = values.magnitude
                units  = values.units
                self.mapping[attr] = pint_pandas.PintArray(magnitude, dtype=units)
            else:
                self.mapping[attr] = [repr(m) for m in values]