#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // For complex number support
#include <pybind11/numpy.h>
#include "single/includes/sphere.cpp"
#include "single/includes/base_class.cpp"


namespace py = pybind11;
using namespace SPHERE;

PYBIND11_MODULE(SphereInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    // Binding for SPHERE::Scatterer class
    py::class_<Scatterer>(module, "SPHERE")
        .def(
            py::init<const double, const complex128, const double, const SOURCE::BaseSource&, size_t>(),
            py::arg("diameter"),
            py::arg("index"),
            py::arg("medium_index"),
            py::arg("source"),
            py::arg("max_order") = 0,
            "Constructor for SPHERE, initializing it with physical and optical properties.")
        .def("get_s1s2", &Scatterer::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters for a sphere.")
        .def("get_fields", &Scatterer::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), py::return_value_policy::move, "Returns the unstructured electromagnetic fields around the sphere.")
        .def("get_full_fields", &Scatterer::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields around the sphere.")
        // Note: Downward are the scattering coeffcients
        .def("an", &Scatterer::get_an_py, py::arg("max_order") = 0, "Returns the an coefficient.")
        .def("bn", &Scatterer::get_bn_py, py::arg("max_order") = 0, "Returns the bn coefficient.")
        .def("cn", &Scatterer::get_cn_py, py::arg("max_order") = 0, "Returns the cn coefficient.")
        .def("dn", &Scatterer::get_dn_py, py::arg("max_order") = 0, "Returns the dn coefficient.")
        // Note: Downward are the efficiencies
        .def_property_readonly("Qsca", &Scatterer::get_Qsca, "Scattering efficiency of the sphere.")
        .def_property_readonly("Qext", &Scatterer::get_Qext, "Extinction efficiency of the sphere.")
        .def_property_readonly("Qabs", &Scatterer::get_Qabs, "Absorption efficiency of the sphere.")
        .def_property_readonly("Qback", &Scatterer::get_Qback, "Backscattering efficiency of the sphere.")
        .def_property_readonly("Qforward", &Scatterer::get_Qforward, "Forward-scattering efficiency of the sphere.")
        .def_property_readonly("Qratio", &Scatterer::get_Qratio, "Ratio of the forward and backward scattering efficiency of the sphere.")
        .def_property_readonly("Qpr", &Scatterer::get_Qpr, "Radiation pressure efficiency of the sphere.")
        // Note: Downward are the cross-sections
        .def_property_readonly("Csca", &Scatterer::get_Csca, "Scattering cross-section of the sphere.")
        .def_property_readonly("Cext", &Scatterer::get_Cext, "Extinction cross-section of the sphere.")
        .def_property_readonly("Cabs", &Scatterer::get_Cabs, "Absorption cross-section of the sphere.")
        .def_property_readonly("Cback", &Scatterer::get_Cback, "Backscattering cross-section of the sphere.")
        .def_property_readonly("Cforward", &Scatterer::get_Cforward, "Forward-scattering efficiency of the sphere.")
        .def_property_readonly("Cratio", &Scatterer::get_Cratio, "Ratio of the forward and backward scattering efficiency of the sphere.")
        .def_property_readonly("Cpr", &Scatterer::get_Cpr, "Radiation pressure cross-section of the sphere.")
        // Note: Downward are the extra parameters
        .def_property_readonly("g", &Scatterer::get_g, "Asymmetry parameter of the sphere.")
        .def_readwrite("area", &Scatterer::area, "Physical cross-sectional area of the sphere.")
        .def_readwrite("size_parameter", &Scatterer::size_parameter, "Size parameter of the sphere scatterer.");
}
