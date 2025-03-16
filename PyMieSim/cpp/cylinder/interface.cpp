#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // Ensure complex numbers are supported
#include "cylinder/cylinder.cpp"

namespace py = pybind11;


PYBIND11_MODULE(interface_cylinder, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<Cylinder>(module, "CYLINDER")
        .def(
            py::init<double, complex128, double, SOURCE::BaseSource&>(),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            "Constructor for CYLINDER, initializing it with physical and optical properties.")

        .def("get_s1s2", &Cylinder::get_s1s2_py, py::arg("phi"), "Calculates and returns the S1 and S2 scattering parameters for a cylinder.")
        .def("get_fields", &Cylinder::get_unstructured_fields_py, py::arg("phi"), py::arg("theta"), py::arg("r"), "Returns the unstructured electromagnetic fields around the cylinder.")
        .def("get_full_fields", &Cylinder::get_full_structured_fields_py, py::arg("sampling"), py::arg("r"), "Returns the full structured electromagnetic fields around the sphere.")
        // Note: Downward are the scattering coeffcients
        .def("a1n", &Cylinder::get_a1n_py, py::arg("max_order") = 0, "Returns the a1n coefficient.")
        .def("b1n", &Cylinder::get_b1n_py, py::arg("max_order") = 0, "Returns the b1n coefficient.")
        .def("a2n", &Cylinder::get_a2n_py, py::arg("max_order") = 0, "Returns the a2n coefficient.")
        .def("b2n", &Cylinder::get_b2n_py, py::arg("max_order") = 0, "Returns the b2n coefficient.")
        // Note: Downward are the efficiencies
        .def_property_readonly("Qsca", &Cylinder::get_Qsca, "Scattering efficiency of the cylinder.")
        .def_property_readonly("Qext", &Cylinder::get_Qext, "Extinction efficiency of the cylinder.")
        .def_property_readonly("Qabs", &Cylinder::get_Qabs, "Absorption efficiency of the cylinder.")
        // Note: Downward are the cross-sections
        .def_property_readonly("Csca", &Cylinder::get_Csca, "Scattering cross-section of the cylinder.")
        .def_property_readonly("Cext", &Cylinder::get_Cext, "Extinction cross-section of the cylinder.")
        .def_property_readonly("Cabs", &Cylinder::get_Cabs, "Absorption cross-section of the cylinder.")
        // Note: Downward are the extra parameters
        .def_property_readonly("g", &Cylinder::get_g, "Asymmetry parameter of the cylinder.")
        .def_readwrite("area", &Cylinder::area, "Physical cross-sectional area of the cylinder.")
        .def_readwrite("size_parameter", &Cylinder::size_parameter, "Size parameter of the cylinder scatterer.");
}
