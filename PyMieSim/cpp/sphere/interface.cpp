#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // For complex number support
#include <pybind11/numpy.h>
#include "sphere/sphere.h"


namespace py = pybind11;

PYBIND11_MODULE(interface_sphere, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    // Binding for Sphere class
    py::class_<Sphere>(module, "SPHERE")
        .def(
            py::init<const double, const complex128, const double, const BaseSource&, size_t>(),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            py::arg("max_order") = 0,
            "Constructor for SPHERE, initializing it with physical and optical properties.")
        .def("get_s1s2", &Sphere::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters for a sphere.")
        .def("get_fields", &Sphere::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), py::return_value_policy::move, "Returns the unstructured electromagnetic fields around the sphere.")
        .def("get_full_fields", &Sphere::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields around the sphere.")
        // Note: Downward are the scattering coeffcients
        .def("an", &Sphere::get_an_py, py::arg("max_order") = 0, "Returns the an coefficient.")
        .def("bn", &Sphere::get_bn_py, py::arg("max_order") = 0, "Returns the bn coefficient.")
        .def("cn", &Sphere::get_cn_py, py::arg("max_order") = 0, "Returns the cn coefficient.")
        .def("dn", &Sphere::get_dn_py, py::arg("max_order") = 0, "Returns the dn coefficient.")
        // Note: Downward are the efficiencies
        .def_property_readonly("Qsca", &Sphere::get_Qsca, "Scattering efficiency of the sphere.")
        .def_property_readonly("Qext", &Sphere::get_Qext, "Extinction efficiency of the sphere.")
        .def_property_readonly("Qabs", &Sphere::get_Qabs, "Absorption efficiency of the sphere.")
        .def_property_readonly("Qback", &Sphere::get_Qback, "Backscattering efficiency of the sphere.")
        .def_property_readonly("Qforward", &Sphere::get_Qforward, "Forward-scattering efficiency of the sphere.")
        .def_property_readonly("Qratio", &Sphere::get_Qratio, "Ratio of the forward and backward scattering efficiency of the sphere.")
        .def_property_readonly("Qpr", &Sphere::get_Qpr, "Radiation pressure efficiency of the sphere.")
        // Note: Downward are the cross-sections
        .def_property_readonly("Csca", &Sphere::get_Csca, "Scattering cross-section of the sphere.")
        .def_property_readonly("Cext", &Sphere::get_Cext, "Extinction cross-section of the sphere.")
        .def_property_readonly("Cabs", &Sphere::get_Cabs, "Absorption cross-section of the sphere.")
        .def_property_readonly("Cback", &Sphere::get_Cback, "Backscattering cross-section of the sphere.")
        .def_property_readonly("Cforward", &Sphere::get_Cforward, "Forward-scattering efficiency of the sphere.")
        .def_property_readonly("Cratio", &Sphere::get_Cratio, "Ratio of the forward and backward scattering efficiency of the sphere.")
        .def_property_readonly("Cpr", &Sphere::get_Cpr, "Radiation pressure cross-section of the sphere.")
        // Note: Downward are the extra parameters
        .def_property_readonly("g", &Sphere::get_g, "Asymmetry parameter of the sphere.")
        .def_readwrite("area", &Sphere::area, "Physical cross-sectional area of the sphere.")
        .def_readwrite("size_parameter", &Sphere::size_parameter, "Size parameter of the sphere scatterer.");
}
