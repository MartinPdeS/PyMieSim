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
        Constructs a table of the scatterer's properties formatted for data visualization.
        This method populates the `mapping` dictionary with user-friendly descriptions and formats of the scatterer properties.

        Returns:
            list: A list of visual representations for each property in the `mapping` dictionary that has been populated.
        """
        for attr in self.attributes:
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
