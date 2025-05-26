#include <pybind11/pybind11.h>
#include "mode_field/mode_field.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_mode_field, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    py::class_<ModeID>(module, "MODEID")
        .def(py::init<>())
        .def_readwrite("mode_family", &ModeID::mode_family, "Mode family (e.g., 'LP', 'HG', 'LG', 'NC')")
        .def_readwrite("number_0", &ModeID::number_0, "First mode number (e.g., azimuthal or x-number)")
        .def_readwrite("number_1", &ModeID::number_1, "Second mode number (e.g., radial or y-number)")
        .def("__repr__", [](const ModeID &m) {
            return "<ModeID: " + m.mode_family + ", " + std::to_string(m.number_0) + ", " + std::to_string(m.number_1) + ">";
        })
        .def("__str__", [](const ModeID &m) {
            return m.mode_family + std::to_string(m.number_0) + std::to_string(m.number_1);
        })
    ;

    py::class_<ModeField>(module, "MODEFIELD")
        .def(py::init<const ModeID &>(), py::arg("mode_id"), "Initialize ModeField with a ModeID")
        .def("get_unstructured", &ModeField::get_unstructured, py::arg("x_coords"), py::arg("y_coords"), "Get unstructured mode field for given coordinates")
        .def("get_structured", &ModeField::get_structured_py, py::arg("x_coords"), py::arg("y_coords"), "Get structured mode field as a NumPy array for given coordinates")
    ;

}