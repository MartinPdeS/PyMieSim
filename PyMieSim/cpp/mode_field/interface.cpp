#include <pybind11/pybind11.h>
#include "mode_field/mode_field.cpp"

namespace py = pybind11;

PYBIND11_MODULE(interface_mode_field, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    module.def(
        "get_LP",
        &get_LP_mode_field_py,
        py::arg("x"),
        py::arg("y"),
        py::arg("azimuthal_number"),
        py::arg("radial_number")
    );


    module.def(
        "get_LG",
        &get_LG_mode_field_py,
        py::arg("x"),
        py::arg("y"),
        py::arg("azimuthal_number"),
        py::arg("radial_number")
    );


    module.def(
        "get_HG",
        &get_HG_mode_field_py,
        py::arg("x"),
        py::arg("y"),
        py::arg("x_number"),
        py::arg("y_number")
    );
}