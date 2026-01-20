#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy
import pint_pandas
from PyOptik.material.base_class import BaseMaterial

from PyMieSim.units import RefractiveIndex
from PyMieSim.binary.interface_experiment import ScattererProperties, MediumProperties


class BaseScatterer:
    """
    Base class for scatterer objects.  This class handles the initialization and setup of
    scatterer parameters for use in PyMieSim simulations.

    """

    def _add_refractive_index(
        self, name: str, refractive_index: RefractiveIndex | BaseMaterial
    ) -> None:
        """
        Determines whether the provided refractive_index is a refractive index (Quantity) or a material (BaseMaterial),
        and returns the corresponding values.

        Parameters
        ----------
        refractive_index : Quantity or BaseMaterial
            The core refractive_index to be assigned, which can either be a refractive index (Quantity) or a material (BaseMaterial).

        Raises
        ------
        ValueError:
            If the provided refractive_index is neither a Quantity (refractive index) nor a BaseMaterial.
        """
        CPPClass = MediumProperties if name == "medium" else ScattererProperties

        if all(isinstance(item, RefractiveIndex) for item in refractive_index):
            return CPPClass(properties=refractive_index.magnitude)

        elif all(isinstance(item, BaseMaterial) for item in refractive_index):
            eval_index = numpy.asarray(
                [m.compute_refractive_index(self.source.wavelength) for m in refractive_index]
            )
            return CPPClass(properties=eval_index)

        raise TypeError(
            "All elements in the list must be of type 'Quantity' or 'BaseMaterial', not a mix of both."
        )

    def _generate_mapping(self) -> None:
        """
        Constructs a table of the scatterer's refractive_index formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer refractive_index.

        Returns
        -------
        list
            A list of visual representations for each refractive_index in the `mapping` dictionary that has been populated.
        """
        for attr in self.attributes:
            values = getattr(self, attr)
            if values is None:
                continue

            if hasattr(values, "magnitude"):
                self.mapping["scatterer:" + attr] = pint_pandas.PintArray(
                    values.magnitude, dtype=values.units
                )
            else:
                self.mapping["scatterer:" + attr] = [repr(m) for m in values]