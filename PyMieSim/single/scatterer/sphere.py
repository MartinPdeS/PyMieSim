#!/usr/bin/env python
# -*- coding: utf-8 -*-


from PyMieSim.single.scatterer.base import BaseScatterer
from PyMieSim.binary.interface_single import SPHERE


class Sphere(SPHERE, BaseScatterer):
    """
    Class representing a homogeneous spherical scatterer.
    """
