#!/usr/bin/env python
# -*- coding: utf-8 -*-

from PyMieSim.mesh import FibonacciMesh


def test_fibonacci_mesh():
    mesh = FibonacciMesh(
        max_angle=1.3,
        sampling=100
    )

    mesh.plot()
