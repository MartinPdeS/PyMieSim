#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support

#include "experiment/includes/scatterer_properties.cpp"

namespace py = pybind11;
typedef std::complex<double> complex128;


PYBIND11_MODULE(ScattererPropertiesInterface, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    py::class_<ScattererProperties>(module, "CppScattererProperties")
        .def(py::init<std::vector<complex128>>(), py::arg("index_properties"))
        .def(py::init<std::vector<std::vector<complex128>>>(), py::arg("material_properties"))
        ;

    py::class_<MediumProperties>(module, "CppMediumProperties")
        .def(py::init<std::vector<double>>(), py::arg("index_properties"))
        .def(py::init<std::vector<std::vector<double>>>(), py::arg("material_properties"))
        ;
}
