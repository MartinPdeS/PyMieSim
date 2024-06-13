#include <pybind11/pybind11.h>
#include "sources.cpp"  // Assuming this includes the necessary DETECTOR definitions

namespace py = pybind11;

PYBIND11_MODULE(SourceInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<SOURCE::BaseSource>(module, "BindedBaseSource");

    py::class_<SOURCE::Planewave, SOURCE::BaseSource>(module, "BindedPlanewave")
        .def(
            py::init<double, std::vector<complex128>, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            "Constructs a Planewave source with specified optical properties. ")

        .def_readonly("wavelength", &SOURCE::Planewave::wavelength, "Wavelength of the source.")
        .def_readonly("jones_vector", &SOURCE::Planewave::jones_vector, "Jones vector of the source.")
        .def_readonly("amplitude", &SOURCE::Planewave::amplitude, "Electric field amplitude of the source.")
        ;

    py::class_<SOURCE::Gaussian, SOURCE::BaseSource>(module, "BindedGaussian")
        .def(
            py::init<double, std::vector<complex128>, double, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            "Constructs a Planewave source with specified optical properties. ")

        .def_readonly("wavelength", &SOURCE::Planewave::wavelength, "Wavelength of the source.")
        .def_readonly("jones_vector", &SOURCE::Planewave::jones_vector, "Jones vector of the source.")
        .def_readonly("amplitude", &SOURCE::Planewave::amplitude, "Electric field amplitude of the source.")
        ;
}
