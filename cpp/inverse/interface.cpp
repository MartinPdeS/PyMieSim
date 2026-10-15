#include "inverse.h"

PYBIND11_MODULE(inverse, module) {
    module.doc() = "Native, dependency-free bounded parameter fitting for PyMieSim models.";

    py::class_<FitParameter>(module, "FitParameter", "A named bounded scalar parameter; values may carry Pint units.")
        .def(py::init<std::string, py::object, py::tuple>(), py::arg("name"), py::arg("initial"), py::arg("bounds"))
        .def_readonly("name", &FitParameter::name)
        .def_readonly("initial", &FitParameter::initial)
        .def_readonly("bounds", &FitParameter::bounds);

    py::class_<MeasuredData>(module, "MeasuredData", "Measured values and optional positive standard uncertainties.")
        .def(py::init<py::object, py::object, std::string>(), py::arg("values"), py::arg("uncertainty") = py::none(),
             py::arg("name") = "observation")
        .def_readonly("values", &MeasuredData::values)
        .def_readonly("uncertainty", &MeasuredData::uncertainty)
        .def_readonly("name", &MeasuredData::name);

    py::class_<ParameterFitResult>(module, "ParameterFitResult",
                                   "Fitted parameters, model predictions, residuals, and convergence diagnostics.")
        .def_readonly("parameters", &ParameterFitResult::parameters)
        .def_readonly("initial_parameters", &ParameterFitResult::initial_parameters)
        .def_readonly("prediction", &ParameterFitResult::prediction)
        .def_readonly("observed", &ParameterFitResult::observed)
        .def_readonly("residuals", &ParameterFitResult::residuals)
        .def_readonly("objective", &ParameterFitResult::objective)
        .def_readonly("success", &ParameterFitResult::success)
        .def_readonly("iterations", &ParameterFitResult::iterations)
        .def_readonly("evaluations", &ParameterFitResult::evaluations)
        .def_readonly("message", &ParameterFitResult::message)
        .def_readonly("parameter_names", &ParameterFitResult::parameter_names)
        .def("summary", &ParameterFitResult::summary)
        .def("__repr__", &ParameterFitResult::repr);

    // Preserve existing imports while exposing more descriptive class names.
    module.attr("Parameter") = module.attr("FitParameter");
    module.attr("Observation") = module.attr("MeasuredData");
    module.attr("FitResult") = module.attr("ParameterFitResult");

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
