#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <string>
#include <vector>

namespace py = pybind11;

class LabeledArray {
public:
    LabeledArray(
        py::object values,
        std::vector<std::string> dims,
        py::dict coords = py::dict(),
        py::dict attrs = py::dict(),
        py::object name = py::none()
    );

    py::array to_numpy() const;
    py::array as_numpy() const;
    py::object as_dataframe() const;
    LabeledArray isel(py::dict indexers) const;
    py::object mean(py::object dimension = py::none()) const;
    py::object plot(
        py::object x = py::none(),
        py::object y = py::none(),
        py::object ax = py::none(),
        py::object cmap = py::str("viridis"),
        py::object std_dimension = py::none(),
        double alpha = 0.4,
        bool show_legend = true,
        bool show = true
    ) const;
    py::tuple dims() const;
    py::tuple shape() const;
    std::string repr() const;

    py::array values_;
    std::vector<std::string> dims_;
    py::dict coords_;
    py::dict attrs_;
    py::object name_;
};
