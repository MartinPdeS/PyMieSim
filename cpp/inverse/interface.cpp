#include "inverse.h"

PYBIND11_MODULE(inverse, module) {
    module.doc() = "Native, dependency-free bounded parameter fitting for PyMieSim models.";

    py::class_<Parameter>(module, "Parameter", "A named bounded scalar parameter; values may carry Pint units.")
        .def(py::init<std::string, py::object, py::tuple>(), py::arg("name"), py::arg("initial"), py::arg("bounds"))
        .def_readonly("name", &Parameter::name)
        .def_readonly("initial", &Parameter::initial)
        .def_readonly("bounds", &Parameter::bounds);

    py::class_<Observation>(module, "Observation", "Measured values and optional positive standard uncertainties.")
        .def(py::init<py::object, py::object, std::string>(), py::arg("values"), py::arg("uncertainty") = py::none(),
             py::arg("name") = "observation")
        .def_readonly("values", &Observation::values)
        .def_readonly("uncertainty", &Observation::uncertainty)
        .def_readonly("name", &Observation::name);

    py::class_<FitResult>(module, "FitResult", "Result returned by ``fit_parameters``.")
        .def_readonly("parameters", &FitResult::parameters)
        .def_readonly("initial_parameters", &FitResult::initial_parameters)
        .def_readonly("prediction", &FitResult::prediction)
        .def_readonly("observed", &FitResult::observed)
        .def_readonly("residuals", &FitResult::residuals)
        .def_readonly("objective", &FitResult::objective)
        .def_readonly("success", &FitResult::success)
        .def_readonly("iterations", &FitResult::iterations)
        .def_readonly("evaluations", &FitResult::evaluations)
        .def_readonly("message", &FitResult::message)
        .def_readonly("parameter_names", &FitResult::parameter_names)
        .def("summary", &FitResult::summary)
        .def("__repr__", &FitResult::repr);

    module.def("fit_parameters", &fit_parameters, py::arg("model"), py::arg("observation"), py::arg("parameters"),
               py::kw_only(), py::arg("max_iterations") = 200, py::arg("initial_step") = .1,
               py::arg("step_tolerance") = 1e-6, py::arg("show_progress") = false,
               R"pbdoc(
                   Fit bounded parameters with deterministic coordinate search.

                   Set ``show_progress=True`` to display the current iteration,
                   objective value, and normalized step size. Progress output is
                   disabled by default.
               )pbdoc");
}
