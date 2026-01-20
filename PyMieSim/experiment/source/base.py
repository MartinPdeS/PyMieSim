#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pint_pandas
from dataclasses import fields


class BaseSource:
    """
    Base class for light sources in PyMieSim experiments.
    """

    def _generate_mapping(self) -> None:
        """
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns:
            list: A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """
        for attr in [f for f in self.attributes if f != "source"]:
            values = getattr(self, attr)

            if values is None:
                continue

            if hasattr(values, "magnitude"):
                magnitude = values.magnitude
                units = values.units
                self.mapping["source:" + attr] = pint_pandas.PintArray(
                    magnitude, dtype=units
                )
            else:
                self.mapping["source:" + attr] = [repr(m) for m in values]

    def _generate_mapping(self) -> None:
        """
        Updates the internal mapping of the source with current parameter values, allowing for visual representation
        of the source's properties in a tabular format (useful for debugging and visualization).

        Attributes like wavelength, polarization, optical power, and NA are included in this mapping.
        """
        self.mapping = {}

        for attr in self.attributes:
            value = getattr(self, attr)
            string = "source:" + attr
            if hasattr(value, "magnitude"):
                self.mapping[string] = pint_pandas.PintArray(value, dtype=value.units)
            else:
                self.mapping[string] = [repr(v) for v in value]