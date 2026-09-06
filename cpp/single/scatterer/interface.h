#pragma once
#include <pybind11/pybind11.h>

void register_base_scatterer(pybind11::module_& module);
void register_sphere(pybind11::module_& module);
void register_coreshell(pybind11::module_& module);
void register_cylinder(pybind11::module_& module);
