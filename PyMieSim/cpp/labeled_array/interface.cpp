#include "labeled_array.h"

#include <pybind11/stl.h>

#include <stdexcept>

PYBIND11_MODULE(labeled_array, module)
{
    module.doc() = "Small native labeled-array container for experiment results.";
    py::class_<LabeledArray>(module, "LabeledArray")
        .def(
            py::init<py::object, std::vector<std::string>, py::dict, py::dict, py::object>(),
            py::arg("values"), py::arg("dims"), py::arg("coords") = py::dict(),
            py::arg("attrs") = py::dict(), py::arg("name") = py::none()
        )
        .def_property_readonly("values", &LabeledArray::to_numpy)
        .def_property_readonly("dims", &LabeledArray::dims)
        .def_readonly("coords", &LabeledArray::coords_)
        .def_readonly("attrs", &LabeledArray::attrs_)
        .def_readonly("name", &LabeledArray::name_)
        .def_property_readonly("shape", &LabeledArray::shape)
        .def_property_readonly("ndim", [](const LabeledArray& self) { return self.values_.ndim(); })
        .def("to_numpy", &LabeledArray::to_numpy)
        .def("as_numpy", &LabeledArray::as_numpy)
        .def("as_dataframe", &LabeledArray::as_dataframe)
        .def("isel", &LabeledArray::isel, py::arg("indexers"))
        .def(
            "mean",
            [](const LabeledArray& self, py::object dim, py::object x) {
                if (!dim.is_none() && !x.is_none()) {
                    throw std::invalid_argument("provide either dim= or x=, not both");
                }
                return self.mean(x.is_none() ? dim : x);
            },
            py::arg("dim") = py::none(), py::arg("x") = py::none()
        )
        .def(
            "plot",
            &LabeledArray::plot,
            py::arg("x") = py::none(), py::arg("y") = py::none(),
            py::arg("ax") = py::none(), py::arg("cmap") = py::str("viridis"),
            py::arg("std") = py::none(), py::arg("alpha") = 0.4,
            py::arg("legend") = true, py::arg("show") = true
        )
        .def("__repr__", &LabeledArray::repr);
}
