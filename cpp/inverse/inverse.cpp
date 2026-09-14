#include "inverse.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <string>
#include <vector>

namespace py = pybind11;

py::object magnitude(py::object value) {
    return py::hasattr(value, "magnitude") ? value.attr("magnitude") : value;
}

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

void validate_optimizer_settings(int max_iterations, double initial_step, double step_tolerance) {
    if (max_iterations < 1)
        throw py::value_error("max_iterations must be a positive integer");
    if (initial_step <= 0 || initial_step > 1)
        throw py::value_error("initial_step must be in (0, 1]");
    if (step_tolerance <= 0)
        throw py::value_error("step_tolerance must be positive");
}

std::vector<Parameter> collect_parameters(py::iterable supplied) {
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
    return parameters;
}

py::dict parameter_mapping(const std::vector<Parameter> &parameters, const std::vector<double> &magnitudes) {
    py::dict values;
    for (size_t index = 0; index < parameters.size(); ++index) {
        const py::object &original = parameters[index].initial;
        const double original_magnitude = scalar(original, "fit parameters and bounds must be scalar values");
        values[parameters[index].name.c_str()] =
            py::hasattr(original, "units") ? (py::float_(magnitudes[index]) * original) / py::float_(original_magnitude)
                                           : py::float_(magnitudes[index]);
    }
    return values;
}

double normalized_step_size(const std::vector<double> &step, const std::vector<double> &widths) {
    double largest = 0;
    for (size_t index = 0; index < step.size(); ++index)
        largest = std::max(largest, step[index] / widths[index]);
    return largest;
}

void report_progress(int iteration, double objective, double normalized_step) {
    py::object line = py::str("iteration {:4d} | objective {:.8g} | step {:.3g}")
                          .attr("format")(iteration, objective, normalized_step);
    py::print(line, py::arg("end") = "\r", py::arg("flush") = true);
}

py::object fit_parameters(py::function model, const Observation &observation, py::iterable supplied,
                          int max_iterations = 200, double initial_step = .1, double step_tolerance = 1e-6,
                          bool show_progress = false) {
    validate_optimizer_settings(max_iterations, initial_step, step_tolerance);
    std::vector<Parameter> parameters = collect_parameters(supplied);
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
    std::vector<double> widths(step.size());
    for (size_t index = 0; index < widths.size(); ++index)
        widths[index] = upper[index] - lower[index];

    int evaluations = 0;
    py::object prediction, residuals;
    double best = 0;
    auto score = [&](const std::vector<double> &values) {
        py::dict passed = parameter_mapping(parameters, values);
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
        const double largest = normalized_step_size(step, widths);
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
            report_progress(iterations, best, normalized_step_size(step, widths));
        }
    }
    score(current);
    py::dict values = parameter_mapping(parameters, current);
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
