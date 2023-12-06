#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy


def NA_to_angle(NA: float) -> numpy.ndarray:
    if NA <= 1.0:
        return numpy.arcsin(NA)
    if NA >= 1.0:
        return numpy.arcsin(NA - 1) + numpy.pi / 2


def angle_space_to_direct(angle_space: numpy.ndarray, k: float,) -> numpy.ndarray:
    RadSpace = numpy.deg2rad(angle_space)

    fourier_space = numpy.sin(RadSpace) * k / (2 * numpy.pi)

    fourier_unit = numpy.abs(fourier_space[1] - fourier_space[0])

    DirectSpace = numpy.fft.fftshift(numpy.fft.fftfreq(angle_space.shape[0], d=fourier_unit))

    return DirectSpace


def direct_space_to_angle(direct_space: numpy.ndarray, k: float) -> numpy.ndarray:
    direct_unit = numpy.abs(direct_space[1] - direct_space[0])

    fourier_space = numpy.fft.fftshift(numpy.fft.fftfreq(direct_space.shape[0], d=direct_unit))

    angle_space = numpy.arcsin(2 * numpy.pi * fourier_space / k)  # conversion spatial frequency to angular space

    if numpy.isnan(angle_space).any():
        raise Exception("Magnification too large.")

    return angle_space * 180 / numpy.pi


def cartesian_to_spherical(x: numpy.ndarray, y: numpy.ndarray, z: numpy.ndarray) -> tuple:
    x = numpy.asarray(x)
    y = numpy.asarray(y)
    z = numpy.asarray(z)

    r = numpy.sqrt(x**2 + y**2 + z**2)
    phi = numpy.arcsin(z / r)
    theta = numpy.arctan2(y, x)
    return r, phi, theta


def spherical_to_cartesian(phi: numpy.ndarray, theta: numpy.ndarray, r: numpy.ndarray = None) -> tuple:
    phi = numpy.asarray(phi)
    theta = numpy.asarray(theta)
    r = r if r is not None else numpy.ones(phi.shape)

    x = r * numpy.cos(phi) * numpy.cos(theta)
    y = r * numpy.cos(phi) * numpy.sin(theta)
    z = r * numpy.sin(phi)
    return x, y, z


def rotate_on_y(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    x, y, z = spherical_to_cartesian(phi=phi, theta=theta)

    xp = x * numpy.cos(angle) + z * numpy.sin(angle)
    yp = y
    zp = z * numpy.cos(angle) - x * numpy.sin(angle)
    return cartesian_to_spherical(x=xp, y=yp, z=zp)


def rotate_on_z(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    x, y, z = spherical_to_cartesian(phi=phi, theta=theta)

    xp = x * numpy.cos(angle) - y * numpy.sin(angle)
    yp = x * numpy.sin(angle) + y * numpy.cos(angle)
    zp = z
    return cartesian_to_spherical(x=xp, y=yp, z=zp)


def rotate_on_x(phi: numpy.ndarray, theta: numpy.ndarray, angle: float) -> tuple:
    x, y, z = spherical_to_cartesian(phi=phi, theta=theta)
    xp = x
    yp = y * numpy.cos(angle) - z * numpy.sin(angle)
    zp = y * numpy.sin(angle) + z * numpy.cos(angle)
    return cartesian_to_spherical(x=xp, y=yp, z=zp)

# -
