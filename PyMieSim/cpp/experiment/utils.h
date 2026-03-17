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


template<typename SetType>
py::list get_polarizationset_representation(const SetType& objects) {
    py::list values;

    for (size_t i = 0; i < objects.size(); ++i) {
        const auto& object = objects[i];

        if (object.is_empty) {
            values.append(py::str("None"));
        }
        else if (!std::isnan(object.angle)) {
            values.append(py::cast(object.angle));
        }
        else if (
            object.jones_vector.size() == 2 &&
            object.jones_vector[0] == complex128(1.0 / std::sqrt(2.0), 0.0) &&
            object.jones_vector[1] == complex128(0.0, 1.0 / std::sqrt(2.0))
        ) {
            values.append(py::str("LeftCircular"));
        }
        else if (
            object.jones_vector.size() == 2 &&
            object.jones_vector[0] == complex128(1.0 / std::sqrt(2.0), 0.0) &&
            object.jones_vector[1] == complex128(0.0, -1.0 / std::sqrt(2.0))
        ) {
            values.append(py::str("RightCircular"));
        }
        else {
            std::ostringstream stream;
            stream << "("
                   << object.jones_vector[0] << ", "
                   << object.jones_vector[1] << ")";
            values.append(py::str(stream.str()));
        }
    }

    return values;
}