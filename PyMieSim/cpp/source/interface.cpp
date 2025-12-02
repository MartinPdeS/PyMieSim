#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "source.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_source, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<BaseSource>(module, "BindedBaseSource")
        .def(py::init<>());

    py::class_<Planewave, BaseSource>(module, "PLANEWAVE")
        .def(
            py::init<double, std::vector<complex128>, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            "Constructs a Planewave source with specified optical properties.")
        .def_readonly("_cpp_wavelength", &Planewave::wavelength, "Wavelength of the source.")
        .def_readonly("_cpp_jones_vector", &Planewave::jones_vector, "Jones vector of the source.")
        .def_readonly("_cpp_amplitude", &Planewave::amplitude, "Electric field amplitude of the source.");

    py::class_<Gaussian, BaseSource>(module, "GAUSSIAN")
        .def(
            py::init<double, std::vector<complex128>, double, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            "Constructs a Gaussian source with specified optical properties.")
        .def(py::init<>())
        .def_readonly("_cpp_wavelength", &Gaussian::wavelength, "Wavelength of the source.")
        .def_readonly("_cpp_jones_vector", &Gaussian::jones_vector, "Jones vector of the source.")
        .def_readonly("_cpp_amplitude", &Gaussian::amplitude, "Electric field amplitude of the source.")
        .def_readonly("_cpp_NA", &Gaussian::NA, "Numerical Aperture of the source.")
        .def_readonly("_cpp_optical_power", &Gaussian::optical_power, "Optical power of the source.");
}
