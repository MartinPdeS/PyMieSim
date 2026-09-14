#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

class ParticleSizeDistribution {
public:
    ParticleSizeDistribution(py::object diameters, py::object number_weights);
    ParticleSizeDistribution(std::vector<double> diameters, std::vector<double> weights, bool normalize);

    py::object diameters() const;
    py::array_t<double> number_fractions() const;

    static ParticleSizeDistribution monodisperse(py::object diameter);
    static ParticleSizeDistribution uniform(py::object minimum, py::object maximum, size_t sampling);
    static ParticleSizeDistribution lognormal(py::object median, double geometric_std, size_t sampling);
    static ParticleSizeDistribution truncated_normal(py::object mean, py::object standard_deviation, py::object minimum,
                                                     py::object maximum, size_t sampling);
    static ParticleSizeDistribution triangular(py::object minimum, py::object mode, py::object maximum,
                                               size_t sampling);
    static ParticleSizeDistribution mixture(py::iterable distributions, py::object number_weights);

    std::string repr() const;

private:
    std::vector<double> diameters_;
    std::vector<double> weights_;

    static void validate_sampling(size_t sampling);
    static std::pair<double, double> bounds(py::object minimum, py::object maximum);
    static std::pair<std::vector<double>, std::vector<double>> legendre_unit(size_t sampling);
};

std::pair<std::vector<double>, std::vector<double>> legendre(size_t order);
std::pair<std::vector<double>, std::vector<double>> hermite(size_t order);
size_t parse_sampling(py::object value);
