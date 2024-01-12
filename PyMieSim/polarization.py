#!/usr/bin/env python
# -*- coding: utf-8 -*-

from collections.abc import Iterable

import numpy


class JonesVector():
    def __init__(self, jones_vector: Iterable) -> None:
        self.jones_vector = numpy.array(jones_vector).astype(complex)

    def __repr__(self):
        return self.jones_vector.__repr__()

    def __add__(self, other):
        if self.jones_vector.ndim == 1 and other.jones_vector.ndim == 1:
            return JonesVector([self.jones_vector, other.jones_vector])

        if self.jones_vector.ndim == 2 and other.jones_vector.ndim == 1:
            return JonesVector([*self.jones_vector, other.jones_vector])

        if self.jones_vector.ndim == 1 and other.jones_vector.ndim == 2:
            return JonesVector([self.jones_vector, *other.jones_vector])

        if self.jones_vector.ndim == 2 and other.jones_vector.ndim == 2:
            return JonesVector([*self.jones_vector, *other.jones_vector])


class RightCircularPolarization(JonesVector):
    def __init__(self):
        super().__init__([1, 1j])


class LeftCircularPolarization(JonesVector):
    def __init__(self):
        super().__init__([1, -1j])


class LinearPolarization(JonesVector):
    def __init__(self, *angle_list: list):

        angle_list = numpy.asarray(angle_list).astype(float)

        if numpy.nan in angle_list:
            raise ValueError("Unpolarized light source is not implemented yet...")

        self.angle_list = numpy.atleast_1d(angle_list)

        jones_vector = numpy.cos(self.angle_list * numpy.pi / 180), numpy.sin(self.angle_list * numpy.pi / 180)

        super().__init__(jones_vector=jones_vector)
# -
