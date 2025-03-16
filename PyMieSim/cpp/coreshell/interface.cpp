
#include "coreshell/coreshell.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_coreshell, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<CoreShell>(module, "CORESHELL")
        .def(py::init<double, double, std::complex<double>, std::complex<double>, double, SOURCE::BaseSource&>(),
             py::arg("core_diameter"),
             py::arg("shell_thickness"),
             py::arg("core_refractive_index"),
             py::arg("shell_refractive_index"),
             py::arg("medium_refractive_index"),
             py::arg("source"),
             "Constructor for CORESHELL, initializing it with physical and optical properties.")

        .def("get_s1s2", &CoreShell::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters for a core-shell.")
        .def("get_fields", &CoreShell::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), py::return_value_policy::move, "Returns the unstructured electromagnetic fields around the core-shell.")
        .def("get_full_fields", &CoreShell::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields around the core-shell.")
        // Note: Downward are the scattering coeffcients
        .def("an", &CoreShell::get_an_py, py::arg("max_order") = 0, "Returns the an coefficient.")
        .def("bn", &CoreShell::get_bn_py, py::arg("max_order") = 0, "Returns the bn coefficient.")
        .def("cn", &CoreShell::get_cn_py, py::arg("max_order") = 0, "Returns the cn coefficient, if applicable.")
        .def("dn", &CoreShell::get_dn_py, py::arg("max_order") = 0, "Returns the dn coefficient, if applicable.")
        // Note: Downward are the efficiencies
        .def_property_readonly("Qsca", &CoreShell::get_Qsca, "Scattering efficiency of the core-shell.")
        .def_property_readonly("Qext", &CoreShell::get_Qext, "Extinction efficiency of the core-shell.")
        .def_property_readonly("Qabs", &CoreShell::get_Qabs, "Absorption efficiency of the core-shell.")
        .def_property_readonly("Qback", &CoreShell::get_Qback, "Backscattering efficiency of the core-shell.")
        .def_property_readonly("Qforward", &CoreShell::get_Qforward, "Forward-scattering efficiency of the core-shell.")
        .def_property_readonly("Qratio", &CoreShell::get_Qratio, "Ratio of the forward and backward scattering efficiency of the core-shell.")
        .def_property_readonly("Qpr", &CoreShell::get_Qpr, "Radiation pressure efficiency of the core-shell.")
        // Note: Downward are the cross-sections
        .def_property_readonly("Csca", &CoreShell::get_Csca, "Scattering cross-section of the core-shell.")
        .def_property_readonly("Cext", &CoreShell::get_Cext, "Extinction cross-section of the core-shell.")
        .def_property_readonly("Cabs", &CoreShell::get_Cabs, "Absorption cross-section of the core-shell.")
        .def_property_readonly("Cback", &CoreShell::get_Cback, "Backscattering cross-section of the core-shell.")
        .def_property_readonly("Cforward", &CoreShell::get_Cforward, "Forward-scattering efficiency of the core-shell.")
        .def_property_readonly("Cratio", &CoreShell::get_Cratio, "Ratio of the forward and backward scattering efficiency of the core-shell.")
        .def_property_readonly("Cpr", &CoreShell::get_Cpr, "Radiation pressure cross-section of the core-shell.")
        // Note: Downward are the extra parameters
        .def_property_readonly("g", &CoreShell::get_g, "Asymmetry parameter of the core-shell.")
        .def_readwrite("area", &CoreShell::area, "Physical cross-sectional area of the core-shell.")
        .def_readwrite("size_parameter", &CoreShell::size_parameter, "Size parameter of the core-shell scatterer.");
}
