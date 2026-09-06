#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <material/material.h>

namespace py = pybind11;

std::shared_ptr<BaseMaterial> parse_material_object(
    const py::object& material_object, const py::object& ureg);
std::shared_ptr<BaseMedium> parse_medium_object(
    const py::object& medium_object, const py::object& ureg);
