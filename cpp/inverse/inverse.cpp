#include "inverse.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <string>
#include <vector>

namespace py = pybind11;

py::object magnitude(py::object value) { return py::hasattr(value, "magnitude") ? value.attr("magnitude") : value; }

double scalar(py::object value, const char *message) {
    py::module_ np = py::module_::import("numpy");
    py::object array = np.attr("asarray")(magnitude(std::move(value)), py::arg("dtype") = "float64");
    if (array.attr("ndim").cast<int>() != 0)
        throw py::value_error(message);
    return array.attr("item")().cast<double>();
}

Parameter::Parameter(std::string name, py::object initial, py::tuple bounds)
    : name(std::move(name)), initial(std::move(initial)), bounds(std::move(bounds)) {
    if (this->name.empty())
        throw py::value_error("parameter name must not be empty");
    if (this->bounds.size() != 2)
        throw py::value_error("bounds must contain lower and upper values");
    const double x = scalar(this->initial, "fit parameters and bounds must be scalar values");
    const double lower = scalar(this->bounds[0], "fit parameters and bounds must be scalar values");
    const double upper = scalar(this->bounds[1], "fit parameters and bounds must be scalar values");
    if (!std::isfinite(x) || !std::isfinite(lower) || !std::isfinite(upper))
        throw py::value_error("parameter must have finite initial and bounds");
    if (lower >= upper)
        throw py::value_error("parameter requires lower bound < upper bound");
    if (x < lower || x > upper)
        throw py::value_error("parameter initial value must lie within bounds");
}

Observation::Observation(py::object values, py::object uncertainty, std::string name)
    : values(std::move(values)), uncertainty(std::move(uncertainty)), name(std::move(name)) {
    py::module_ np = py::module_::import("numpy");
    py::object array = np.attr("asarray")(magnitude(this->values), py::arg("dtype") = "float64");
    if (array.attr("size").cast<size_t>() == 0 || !np.attr("all")(np.attr("isfinite")(array)).cast<bool>())
        throw py::value_error("observation values must be nonempty and finite");
    if (!this->uncertainty.is_none()) {
        py::object sigma = np.attr("asarray")(magnitude(this->uncertainty), py::arg("dtype") = "float64");
        try {
            np.attr("broadcast_to")(sigma, array.attr("shape"));
        } catch (py::error_already_set &) {
            throw py::value_error("uncertainty must be broadcastable to observation shape");
        }
        if (!np.attr("all")(np.attr("isfinite")(sigma)).cast<bool>() ||
            np.attr("any")(np.attr("less_equal")(sigma, 0)).cast<bool>())
            throw py::value_error("uncertainty must be finite and positive");
    }
}

FitResult::FitResult() : parameter_names(0) {}

std::string FitResult::repr() const {
    return "FitResult(" + std::string(success ? "converged" : "stopped") + ", objective=" + std::to_string(objective) +
           ", parameters=" + py::str(parameters).cast<std::string>() + ")";
}

std::string FitResult::summary() const {
    return message + " (" + std::to_string(evaluations) +
           " model evaluations)\nobjective = " + std::to_string(objective);
}

py::object fit_parameters(py::function model, const Observation &observation, py::iterable supplied,
                          int max_iterations = 200, double initial_step = .1, double step_tolerance = 1e-6,
                          bool show_progress = false) {
    if (max_iterations < 1 || initial_step <= 0 || initial_step > 1 || step_tolerance <= 0)
        throw py::value_error("invalid optimizer settings");
    std::vector<Parameter> parameters;
    py::set names;
    for (py::handle item : supplied) {
        auto parameter = item.cast<Parameter>();
        if (names.contains(parameter.name.c_str()))
            throw py::value_error("parameters must be nonempty and have unique names");
        names.add(parameter.name.c_str());
        parameters.push_back(std::move(parameter));
    }
    if (parameters.empty())
        throw py::value_error("parameters must be nonempty and have unique names");
    py::module_ np = py::module_::import("numpy");
    py::object observed = np.attr("asarray")(magnitude(observation.values), py::arg("dtype") = "float64");
    py::object sigma = observation.uncertainty.is_none()
                           ? np.attr("ones_like")(observed)
                           : np.attr("broadcast_to")(
                                 np.attr("asarray")(magnitude(observation.uncertainty), py::arg("dtype") = "float64"),
                                 observed.attr("shape"));
    std::vector<double> lower, upper, current, step;
    py::dict initial;
    for (const auto &p : parameters) {
        double lo = scalar(p.bounds[0], "fit parameters and bounds must be scalar values"),
               hi = scalar(p.bounds[1], "fit parameters and bounds must be scalar values");
        lower.push_back(lo);
        upper.push_back(hi);
        current.push_back(scalar(p.initial, "fit parameters and bounds must be scalar values"));
        step.push_back((hi - lo) * initial_step);
        initial[p.name.c_str()] = p.initial;
    }
    int evaluations = 0;
    py::object prediction, residuals;
    double best = 0;
    auto score = [&](const std::vector<double> &values) {
        py::dict passed;
        for (size_t i = 0; i < parameters.size(); ++i) {
            auto &original = parameters[i].initial;
            double base = scalar(original, "fit parameters and bounds must be scalar values");
            passed[parameters[i].name.c_str()] = py::hasattr(original, "units")
                                                     ? (py::float_(values[i]) * original) / py::float_(base)
                                                     : py::float_(values[i]);
        }
        prediction = np.attr("asarray")(magnitude(model(passed)), py::arg("dtype") = "float64");
        if (!py::bool_(prediction.attr("shape").equal(observed.attr("shape"))))
            throw py::value_error("model returned shape different from observation");
        residuals = (prediction - observed) / sigma;
        ++evaluations;
        return np.attr("dot")(residuals.attr("ravel")(), residuals.attr("ravel")()).cast<double>();
    };
    best = score(current);
    int iterations = 0;
    while (iterations < max_iterations) {
        double largest = 0;
        for (size_t i = 0; i < step.size(); ++i)
            largest = std::max(largest, step[i] / (upper[i] - lower[i]));
        if (largest < step_tolerance)
            break;
        ++iterations;
        bool improved = false;
        for (size_t i = 0; i < current.size(); ++i)
            for (double direction : {-1., 1.}) {
                auto candidate = current;
                candidate[i] = std::max(lower[i], std::min(upper[i], candidate[i] + direction * step[i]));
                if (candidate[i] == current[i])
                    continue;
                double value = score(candidate);
                if (value < best) {
                    current = candidate;
                    best = value;
                    improved = true;
                }
            }
        if (!improved)
            for (auto &value : step)
                value *= .5;

        if (show_progress) {
            double normalized_step = 0;
            for (size_t i = 0; i < step.size(); ++i)
                normalized_step = std::max(normalized_step, step[i] / (upper[i] - lower[i]));

            py::object line = py::str("iteration {:4d} | objective {:.8g} | step {:.3g}")
                                  .attr("format")(iterations, best, normalized_step);
            py::print(line, py::arg("end") = "\r", py::arg("flush") = true);
        }
    }
    score(current);
    py::dict values;
    for (size_t i = 0; i < parameters.size(); ++i) {
        auto &original = parameters[i].initial;
        double base = scalar(original, "fit parameters and bounds must be scalar values");
        values[parameters[i].name.c_str()] = py::hasattr(original, "units")
                                                 ? (py::float_(current[i]) * original) / py::float_(base)
                                                 : py::float_(current[i]);
    }
    auto result = FitResult{};
    result.parameters = values;
    result.initial_parameters = initial;
    result.prediction = prediction;
    result.observed = observed;
    result.residuals = residuals;
    result.objective = best;
    result.iterations = iterations;
    result.evaluations = evaluations;
    result.success = iterations < max_iterations;
    result.message = result.success ? "converged" : "maximum iterations reached";
    if (show_progress) {
        py::object line =
            py::str("{} after {} iterations | objective {:.8g}").attr("format")(result.message, iterations, best);
        py::print(line);
    }
    py::tuple names_tuple(parameters.size());
    for (size_t i = 0; i < parameters.size(); ++i)
        names_tuple[i] = parameters[i].name;
    result.parameter_names = names_tuple;
    return py::cast(result);
}
