#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // For complex number support
#include <pybind11/numpy.h>
#include "scatterer/sphere/sphere.h"
#include "scatterer/coreshell/coreshell.h"
#include "scatterer/cylinder/cylinder.h"



namespace py = pybind11;

PYBIND11_MODULE(interface_scatterer, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<BaseScatterer>(module, "BASESCATTERER")
        // Note: Downward are the scattering fields
        .def("get_s1s2", &BaseScatterer::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters for a scatterer.")
        .def("get_fields", &BaseScatterer::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), py::return_value_policy::move, "Returns the unstructured electromagnetic fields around the scatterer.")
        .def("get_full_fields", &BaseScatterer::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields around the scatterer.")
        // Note: Downward are the scattering coeffcients
        .def("an", &BaseScatterer::get_an_list_py, "Returns the an coefficient.")
        .def("bn", &BaseScatterer::get_bn_list_py, "Returns the bn coefficient.")
        .def("cn", &BaseScatterer::get_cn_list_py, "Returns the cn coefficient.")
        .def("dn", &BaseScatterer::get_dn_list_py, "Returns the dn coefficient.")
        .def("get_coefficient", &BaseScatterer::get_coefficient_py, py::arg("type"), py::arg("order"), "Returns the dn coefficient.")
        // Note: Downward are the efficiencies
        .def_property_readonly("Qsca", &BaseScatterer::get_Qsca, "Scattering efficiency of the scatterer.")
        .def_property_readonly("Qext", &BaseScatterer::get_Qext, "Extinction efficiency of the scatterer.")
        .def_property_readonly("Qabs", &BaseScatterer::get_Qabs, "Absorption efficiency of the scatterer.")
        .def_property_readonly("Qback", &BaseScatterer::get_Qback, "Backscattering efficiency of the scatterer.")
        .def_property_readonly("Qforward", &BaseScatterer::get_Qforward, "Forward-scattering efficiency of the scatterer.")
        .def_property_readonly("Qratio", &BaseScatterer::get_Qratio, "Ratio of the forward and backward scattering efficiency of the scatterer.")
        .def_property_readonly("Qpr", &BaseScatterer::get_Qpr, "Radiation pressure efficiency of the scatterer.")
        // Note: Downward are the cross-sections
        .def_property_readonly("Csca", &BaseScatterer::get_Csca, "Scattering cross-section of the scatterer.")
        .def_property_readonly("Cext", &BaseScatterer::get_Cext, "Extinction cross-section of the scatterer.")
        .def_property_readonly("Cabs", &BaseScatterer::get_Cabs, "Absorption cross-section of the scatterer.")
        .def_property_readonly("Cback", &BaseScatterer::get_Cback, "Backscattering cross-section of the scatterer.")
        .def_property_readonly("Cforward", &BaseScatterer::get_Cforward, "Forward-scattering efficiency of the scatterer.")
        .def_property_readonly("Cratio", &BaseScatterer::get_Cratio, "Ratio of the forward and backward scattering efficiency of the scatterer.")
        .def_property_readonly("Cpr", &BaseScatterer::get_Cpr, "Radiation pressure cross-section of the scatterer.")
        // Note: Downward are the extra parameters
        .def_property_readonly("g", &BaseScatterer::get_g, "Asymmetry parameter of the scatterer.")
        .def_readwrite("area", &BaseScatterer::area, "Physical cross-sectional area of the scatterer.")
        .def_readwrite("size_parameter", &BaseScatterer::size_parameter, "Size parameter of the scatterer.")
        ;


    // Binding for Sphere class
    py::class_<Sphere, BaseScatterer>(module, "SPHERE")
        .def(
            py::init<const double, const complex128, const double, const BaseSource&, size_t>(),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            py::arg("max_order") = 0,
            "Constructor for SPHERE, initializing it with physical and optical properties.")
        ;

    // Binding for CoreShell class
    py::class_<CoreShell, BaseScatterer>(module, "CORESHELL")
        .def(py::init<double, double, std::complex<double>, std::complex<double>, double, BaseSource&>(),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            "Constructor for CORESHELL, initializing it with physical and optical properties.")
        ;

    // Binding for Cylinder class
    py::class_<Cylinder, BaseScatterer>(module, "CYLINDER")
        .def(
            py::init<double, complex128, double, BaseSource&>(),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            "Constructor for CYLINDER, initializing it with physical and optical properties.")
        // Note: Downward are the scattering coeffcients
        .def("a1n", &Cylinder::get_a1n_list_py, "Returns the a1n coefficient.")
        .def("b1n", &Cylinder::get_b1n_list_py, "Returns the b1n coefficient.")
        .def("a2n", &Cylinder::get_a2n_list_py, "Returns the a2n coefficient.")
        .def("b2n", &Cylinder::get_b2n_list_py, "Returns the b2n coefficient.")
        ;


}

