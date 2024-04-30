#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // For complex number support
#include <pybind11/numpy.h>
#include "sphere.cpp"


namespace py = pybind11;
using namespace SPHERE;

PYBIND11_MODULE(SphereInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    // Binding for SPHERE::Scatterer class
    py::class_<Scatterer>(module, "SPHERE")
        .def(py::init<double, double, double, std::complex<double>, double, std::vector<complex128>>(),
             py::arg("wavelength"),
             py::arg("amplitude"),
             py::arg("diameter"),
             py::arg("index"),
             py::arg("medium_index"),
             py::arg("jones_vector"),
             "Constructor for SPHERE, initializing it with physical and optical properties.")

        .def(py::init<>(), "Default constructor for SPHERE scatterer.")
        .def("get_s1s2", &Scatterer::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters for a sphere.")
        .def("get_fields", &Scatterer::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), py::return_value_policy::move, "Returns the unstructured electromagnetic fields around the sphere.")
        .def("get_full_fields", &Scatterer::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields around the sphere.")
        .def("an", py::overload_cast<>(&Scatterer::get_an_py), "Returns the an coefficient.")
        .def("bn", py::overload_cast<>(&Scatterer::get_bn_py), "Returns the bn coefficient.")
        .def("cn", py::overload_cast<>(&Scatterer::get_cn_py), "Returns the cn coefficient, if applicable.")
        .def("dn", py::overload_cast<>(&Scatterer::get_dn_py), "Returns the dn coefficient, if applicable.")
        .def_property_readonly("Qsca", &Scatterer::get_Qsca, "Scattering efficiency of the sphere.")
        .def_property_readonly("Qext", &Scatterer::get_Qext, "Extinction efficiency of the sphere.")
        .def_property_readonly("Qabs", &Scatterer::get_Qabs, "Absorption efficiency of the sphere.")
        .def_property_readonly("Qback", &Scatterer::get_Qback, "Backscattering efficiency of the sphere.")
        // Note: Qforward is commented out; assuming it's either not implemented or was an oversight.
        .def_property_readonly("Qpr", &Scatterer::get_Qpr, "Radiation pressure efficiency of the sphere.")
        .def_property_readonly("Csca", &Scatterer::get_Csca, "Scattering cross-section of the sphere.")
        .def_property_readonly("Cext", &Scatterer::get_Cext, "Extinction cross-section of the sphere.")
        .def_property_readonly("Cabs", &Scatterer::get_Cabs, "Absorption cross-section of the sphere.")
        .def_property_readonly("Cback", &Scatterer::get_Cback, "Backscattering cross-section of the sphere.")
        .def_property_readonly("Cpr", &Scatterer::get_Cpr, "Radiation pressure cross-section of the sphere.")
        .def_property_readonly("g", &Scatterer::get_g, "Asymmetry parameter of the sphere.")
        .def_readwrite("area", &Scatterer::area, "Physical cross-sectional area of the sphere.")
        .def_readwrite("size_parameter", &Scatterer::size_parameter, "Size parameter of the sphere scatterer.");
}
