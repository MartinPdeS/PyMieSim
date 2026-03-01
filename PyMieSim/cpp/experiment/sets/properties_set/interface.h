#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support

#include <pint/pint.h>
#include "./properties_set.h"

namespace py = pybind11;
typedef std::complex<double> complex128;



void register_properties_set(py::module& module) {
    py::object ureg = get_shared_ureg();

    py::class_<ScattererProperties, std::shared_ptr<ScattererProperties>>(module, "ScattererProperties")
        .def(
            py::init<std::vector<complex128>>(),
            py::arg("properties")
        )
        .def(
            py::init<std::vector<std::vector<complex128>>>(),
            py::arg("properties")
        )
        ;

    py::class_<MediumProperties, std::shared_ptr<MediumProperties>>(module, "MediumProperties")
        .def(
            py::init<std::vector<double>>(),
            py::arg("properties")
        )
        .def(
            py::init<std::vector<std::vector<double>>>(),
            py::arg("properties")
        )
        ;
}

// -
