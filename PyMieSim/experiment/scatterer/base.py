#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyOptik.base_class import BaseMaterial
from pydantic import field_validator
from typing import Optional, List
from pydantic.dataclasses import dataclass
from dataclasses import fields
from pydantic import ConfigDict
from PyMieSim.units import Quantity
import pint_pandas

import numpy
from PyMieSim.experiment.source.base import BaseSource

config_dict = ConfigDict(
    kw_only=True,
    slots=True,
    extra='forbid',
    arbitrary_types_allowed=True
)


@dataclass(config=config_dict, kw_only=True)
class BaseScatterer():
    """
    Base class for scatterer objects. This class handles the initialization and setup of
    scatterer parameters for use in PyMieSim simulations.

    """
    source: BaseSource
    medium_property: List[BaseMaterial] | List[Quantity]


    mapping = None
    binding_kwargs = None
    binding = None

    @field_validator('property', 'medium_property', mode='plain')
    def _validate_properties(cls, value):
        """Ensure that arrays are properly converted to numpy arrays."""

        return numpy.atleast_1d(value)

    def _assign_index_or_material(self, element: str, property: Quantity | BaseMaterial) -> tuple[Quantity | None, BaseMaterial | None]:
        """
        Determines whether the provided property is a refractive index (Quantity) or a material (BaseMaterial),
        and returns the corresponding values.

        Parameters:
        ----------
        property : Quantity or BaseMaterial
            The core property to be assigned, which can either be a refractive index (Quantity) or a material (BaseMaterial).

        Returns:
        -------
        tuple[Quantity | None, BaseMaterial | None]
            A tuple where the first element is the refractive index (Quantity) if provided, otherwise None.
            The second element is the material (BaseMaterial) if provided, otherwise None.

        Raises:
        ------
        ValueError:
            If the provided property is neither a Quantity (refractive index) nor a BaseMaterial.
        """
        if all(isinstance(item, Quantity) for item in property):
            self.binding_kwargs[f'{element}index'] = property

        elif all(isinstance(item, BaseMaterial) for item in property):
            eval_index = numpy.asarray([m.compute_refractive_index(self.source.wavelength.to_base_units().magnitude) for m in property])
            self.binding_kwargs[f'{element}material'] = eval_index

        else:
            raise TypeError("All elements in the list must be of type 'Quantity' or 'BaseMaterial', not a mix of both.")

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