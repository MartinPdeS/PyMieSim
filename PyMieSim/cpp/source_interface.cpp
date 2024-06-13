#include <pybind11/pybind11.h>
#include "sources.h"

namespace py = pybind11;
using namespace SOURCE;

PYBIND11_MODULE(SourceInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<BaseSource>(module, "BindedBaseSource");

    py::class_<Planewave>(module, "BindedPlanewave")
        .def(
            py::init<double, std::vector<complex128>, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            "Constructs a Planewave source with specified optical properties. ")

        .def_readonly("wavelength", &Planewave::wavelength, "Wavelength of the source.")
        .def_readonly("jones_vector", &Planewave::jones_vector, "Jones vector of the source.")
        .def_readonly("amplitude", &Planewave::amplitude, "Electric field amplitude of the source.")
        ;

    py::class_<Gaussian, BaseSource>(module, "BindedGaussian")
        .def(
            py::init<double, std::vector<complex128>, double, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            "Constructs a Gaussian source with specified optical properties. ")

        .def(py::init<>())

        .def_readonly("wavelength", &Gaussian::wavelength, "Wavelength of the source.")
        .def_readonly("jones_vector", &Gaussian::jones_vector, "Jones vector of the source.")
        .def_readonly("amplitude", &Gaussian::amplitude, "Electric field amplitude of the source.")
        ;
}
