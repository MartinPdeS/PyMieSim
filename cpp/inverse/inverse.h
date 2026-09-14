#pragma once

#include <pybind11/pybind11.h>
#include <string>

namespace py = pybind11;

class Parameter {
  public:
    Parameter(std::string name, py::object initial, py::tuple bounds);
    std::string name;
    py::object initial;
    py::tuple bounds;
};

class Observation {
  public:
    Observation(py::object values, py::object uncertainty, std::string name);
    py::object values;
    py::object uncertainty;
    std::string name;
};

class FitResult {
  public:
    FitResult();
    std::string repr() const;
    std::string summary() const;
    py::dict parameters;
    py::dict initial_parameters;
    py::object prediction;
    py::object observed;
    py::object residuals;
    double objective;
    bool success;
    int iterations;
    int evaluations;
    std::string message;
    py::tuple parameter_names;
};

py::object fit_parameters(py::function model, const Observation &observation, py::iterable parameters,
                          int max_iterations, double initial_step, double step_tolerance);
