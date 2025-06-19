#include <pybind11/pybind11.h>
#include "mode_field.h"


void register_mode_field(py::module_& module) {
    // ------------------ Bindings for ModeField ------------------
    pybind11::class_<ModeField>(module, "MODEFIELD")
        .def("_cpp_get_unstructured",
            &ModeField::get_unstructured,
            pybind11::arg("x_coords"),
            pybind11::arg("y_coords"),
            pybind11::return_value_policy::move,
            R"pbdoc(
                Generate a structured scalar field as a numpy array.

                Parameters
                ----------
                sampling : int
                    The sampling rate for the scalar field. Default is 100.

                Returns
                -------
                numpy.ndarray
                    A 2D array representing the structured scalar field.
            )pbdoc"
        )
    ;
}

