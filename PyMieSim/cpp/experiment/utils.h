#pragma once

#include <memory>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <experiment/material_set/material_set.h>


template<typename SetType>
py::list get_materialset_representation(const SetType& objects) {
    py::list values;

    for (size_t i = 0; i < objects.size(); ++i) {
        py::object object = py::cast(objects[i]);

        if (py::hasattr(object, "name")) {
            values.append(object.attr("name"));
        }
        else {
            values.append(py::cast(objects[i]->refractive_index));
        }
    }

    return values;
}