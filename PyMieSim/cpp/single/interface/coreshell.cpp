#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // Ensure complex numbers are supported
#include "single/includes/coreshell.cpp"
#include <single/includes/base_class.cpp>

namespace py = pybind11;
using namespace CORESHELL;

PYBIND11_MODULE(CoreShellInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<Scatterer>(module, "CORESHELL")
        .def(py::init<double, double, std::complex<double>, std::complex<double>, double, SOURCE::BaseSource&>(),
             py::arg("core_diameter"),
             py::arg("shell_thickness"),
             py::arg("core_index"),
             py::arg("shell_index"),
             py::arg("medium_index"),
             py::arg("source"),
             "Constructor for CORESHELL, initializing it with physical and optical properties.")

        .def("get_s1s2", &Scatterer::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters for a core-shell.")
        .def("get_fields", &Scatterer::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), py::return_value_policy::move, "Returns the unstructured electromagnetic fields around the core-shell.")
        .def("get_full_fields", &Scatterer::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields around the core-shell.")
        // Note: Downward are the scattering coeffcients
        .def("an", &Scatterer::get_an_py, py::arg("max_order") = 0, "Returns the an coefficient.")
        .def("bn", &Scatterer::get_bn_py, py::arg("max_order") = 0, "Returns the bn coefficient.")
        .def("cn", &Scatterer::get_cn_py, py::arg("max_order") = 0, "Returns the cn coefficient, if applicable.")
        .def("dn", &Scatterer::get_dn_py, py::arg("max_order") = 0, "Returns the dn coefficient, if applicable.")
        // Note: Downward are the efficiencies
        .def_property_readonly("Qsca", &Scatterer::get_Qsca, "Scattering efficiency of the core-shell.")
        .def_property_readonly("Qext", &Scatterer::get_Qext, "Extinction efficiency of the core-shell.")
        .def_property_readonly("Qabs", &Scatterer::get_Qabs, "Absorption efficiency of the core-shell.")
        .def_property_readonly("Qback", &Scatterer::get_Qback, "Backscattering efficiency of the core-shell.")
        .def_property_readonly("Qforward", &Scatterer::get_Qforward, "Forward-scattering efficiency of the core-shell.")
        .def_property_readonly("Qratio", &Scatterer::get_Qratio, "Ratio of the forward and backward scattering efficiency of the core-shell.")
        .def_property_readonly("Qpr", &Scatterer::get_Qpr, "Radiation pressure efficiency of the core-shell.")
        // Note: Downward are the cross-sections
        .def_property_readonly("Csca", &Scatterer::get_Csca, "Scattering cross-section of the core-shell.")
        .def_property_readonly("Cext", &Scatterer::get_Cext, "Extinction cross-section of the core-shell.")
        .def_property_readonly("Cabs", &Scatterer::get_Cabs, "Absorption cross-section of the core-shell.")
        .def_property_readonly("Cback", &Scatterer::get_Cback, "Backscattering cross-section of the core-shell.")
        .def_property_readonly("Cforward", &Scatterer::get_Cforward, "Forward-scattering efficiency of the core-shell.")
        .def_property_readonly("Cratio", &Scatterer::get_Cratio, "Ratio of the forward and backward scattering efficiency of the core-shell.")
        .def_property_readonly("Cpr", &Scatterer::get_Cpr, "Radiation pressure cross-section of the core-shell.")
        // Note: Downward are the extra parameters
        .def_property_readonly("g", &Scatterer::get_g, "Asymmetry parameter of the core-shell.")
        .def_readwrite("area", &Scatterer::area, "Physical cross-sectional area of the core-shell.")
        .def_readwrite("size_parameter", &Scatterer::size_parameter, "Size parameter of the core-shell scatterer.");
}
