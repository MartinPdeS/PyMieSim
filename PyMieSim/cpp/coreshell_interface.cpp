#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // Ensure complex numbers are supported
#include "core_shell.cpp"

namespace py = pybind11;
using namespace CORESHELL;

PYBIND11_MODULE(CoreShellInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    // Binding for CoreShell::State class
    py::class_<State>(module, "CppCoreShellState")
        .def(py::init<>()); 

    // Binding for CoreShell::Scatterer class
    py::class_<Scatterer>(module, "CORESHELL")
        .def(py::init<double, double, double, double, std::complex<double>, std::complex<double>, double, CVector>(),
             py::arg("wavelength"),
             py::arg("amplitude"),
             py::arg("core_diameter"),
             py::arg("shell_width"),
             py::arg("core_index"),
             py::arg("shell_index"),
             py::arg("n_medium"),
             py::arg("jones_vector"),
             "Constructor for CORESHELL, initializing it with physical and optical properties.")

        .def("get_s1s2", &Scatterer::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters.")
        .def("get_fields", &Scatterer::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), "Returns the unstructured electromagnetic fields.")
        .def("get_full_fields", &Scatterer::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields.")
        .def("an", py::overload_cast<>(&Scatterer::get_an_py), "Returns the an coefficient.")
        .def("bn", py::overload_cast<>(&Scatterer::get_bn_py), "Returns the bn coefficient.")
        .def_property_readonly("Qsca", &Scatterer::get_Qsca, "Scattering efficiency.")
        .def_property_readonly("Qext", &Scatterer::get_Qext, "Extinction efficiency.")
        .def_property_readonly("Qabs", &Scatterer::get_Qabs, "Absorption efficiency.")
        .def_property_readonly("Qback", &Scatterer::get_Qback, "Backscattering efficiency.")
        .def_property_readonly("Qpr", &Scatterer::get_Qpr, "Radiation pressure efficiency.")
        .def_property_readonly("Csca", &Scatterer::get_Csca, "Scattering cross-section.")
        .def_property_readonly("Cext", &Scatterer::get_Cext, "Extinction cross-section.")
        .def_property_readonly("Cabs", &Scatterer::get_Cabs, "Absorption cross-section.")
        .def_property_readonly("Cback", &Scatterer::get_Cback, "Backscattering cross-section.")
        .def_property_readonly("Cpr", &Scatterer::get_Cpr, "Radiation pressure cross-section.")
        .def_property_readonly("g", &Scatterer::get_g, "Asymmetry parameter.")
        .def_readwrite("state", &Scatterer::state, "State of the scatterer.")
        .def_readwrite("area", &Scatterer::area, "Physical cross-sectional area.")
        .def_readwrite("size_parameter", &Scatterer::size_parameter, "Size parameter of the scatterer.");
}
