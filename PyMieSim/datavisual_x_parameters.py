#!/usr/bin/env python
# -*- coding: utf-8 -*-


# Source
wavelength = dict(
    name='wavelength',
    long_label=r'$\lambda$',
    format=".1e",
    unit="m",
    short_label=r'$\lambda$'
)

amplitude = dict(
    long_label='Ampltiude',
    format=".1e",
    unit="w.m⁻¹",
    short_label=r'E$_{0}$'
)

linear_polarization = dict(
    name='polarization',
    long_label='Linear polarization',
    format=".1f",
    unit="Deg",
    short_label=r'$\hat{P}$'
)

# Detector
scalarfield = dict(
    name='Field',
    long_label='Coupling field',
    format="",
    unit="",
    short_label='C. F.'
)

NA = dict(
    name='numerical_aperture',
    long_label='Numerical aperture',
    format=".3f",
    unit="rad",
    short_label='NA'
)

phi_offset = dict(
    name='phi_offset',
    long_label='Phi angle',
    format="03.1f",
    unit="deg",
    short_label=r'$\phi_{offset}$'
)

gamma_offset = dict(
    name='gamma_offset',
    long_label='Gamma angle',
    format="03.1f",
    unit="deg",
    short_label=r'$\gamma_{offset}$'
)

polarization_filter = dict(
    name='polarization_filter',
    long_label=r'Polarization filter',
    format="03.1f",
    unit="deg",
    short_label=r'f$_{pol}$',
)


# Scatterer

diameter = dict(
    name='diameter',
    format=".2e",
    unit="m",
    long_label='Scatterer diameter',
    short_label='diameter'
)

material = dict(
    name='material',
    format="",
    unit="",
    long_label='Scatterer material',
    short_label='material'
)

index = dict(
    name='index',
    format="",
    unit="1",
    long_label='Refractive index',
    short_label='index'
)

core_diameter = dict(
    name='diameter',
    format=".2e",
    unit="m",
    long_label=r'Core diameter',
    short_label=r'core$_{diameter}$'
)

shell_width = dict(
    name='Shell diameter',
    format=".2e",
    unit="m",
    long_label=r'Shell diameter',
    short_label=r'shell$_{diameter}$'
)

core_material = dict(
    name='core_material',
    format="",
    unit="",
    long_label=r'Core material',
    short_label=r'core$_{material}$'
)

shell_material = dict(
    name='shell_material',
    format="",
    unit="",
    long_label=r'Shell material',
    short_label=r'shell$_{material}$',
)

core_index = dict(
    name='core_index',
    format="",
    unit="1",
    long_label=r'Core index',
    short_label=r'core$_{index}$'
)

shell_index = dict(
    name='shell_index',
    format="",
    unit="1",
    long_label=r'Shell index',
    short_label=r'shell$_{index}$'
)


n_medium = dict(
    name='n_medium',
    long_label=r'Refractive index of medium',
    format=".2f",
    unit="1",
    short_label=r'n$_{medium}$'
)


# -
