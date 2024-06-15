#!/usr/bin/env python
# -*- coding: utf-8 -*-

from DataVisual import units

# Length units
diameter = units.Length(long_label='Scatterer diameter', short_label='diameter', string_format='.1f')
core_diameter = units.Length(long_label='Core diameter', short_label='core_diameter', string_format='.1f')
shell_width = units.Length(long_label='Shell width', short_label='shell_width', string_format='.1f')
wavelength = units.Length(long_label='Wavelength', short_label=r'$\lambda$', string_format='.1f')


# Angle units
polarization = units.Degree(long_label=r'Polarization angle', short_label=r'pol.')
rotation = units.Degree(long_label='Rotation angle', short_label='rot', string_format='.1f')
phi_offset = units.Degree(long_label='Phi angle', short_label=r'phi', use_prefix=False, string_format='.1f')
gamma_offset = units.Degree(long_label='Gamma angle', short_label=r'gamma', use_prefix=False, string_format='.1f')
polarization_filter = units.Degree(long_label=r'Polarization filter', short_label=r'f$_{pol}$', use_prefix=False, string_format='.1f')

# Index units
NA = units.Index(long_label='Numerical aperture', short_label='NA', use_prefix=False, string_format=".2f")
NA_source = units.Index(long_label='Source numerical aperture', short_label='source NA', use_prefix=False, string_format=".2f")
index = units.Index(long_label='Refractive index', short_label='index', use_prefix=False, string_format=".2f")
core_index = units.Index(long_label='Core refractive index', short_label='core index', use_prefix=False, string_format=".2f")
shell_index = units.Index(long_label='Shell refractive index', short_label='shell index', use_prefix=False, string_format=".2f")
medium_index = units.Index(long_label='Medium refractive index', short_label='medium index', use_prefix=False, string_format=".2f")


# Custom units
material = units.Custom(long_label='Material', short_label='material')
core_material = units.Custom(long_label='Core material', short_label='core material')
shell_material = units.Custom(long_label='Shell material', short_label='shell material')
medium_material = units.Custom(long_label='Medium material', short_label='medium material')
optical_power = units.Custom(long_label='Optical power', short_label='power')

# Custom units
sampling = units.Custom(long_label='Sampling', short_label='sampling')
mode_number = units.Custom(long_label='Mode number', short_label='mode')

# -
