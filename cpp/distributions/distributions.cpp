#include "distributions.h"
#include <pint/pint.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace py = pybind11;

std::pair<std::vector<double>, std::vector<double>> legendre(size_t order) {
    if (order == 0)
        throw std::invalid_argument("quadrature order must be positive");
    std::vector<double> nodes(order), weights(order);
    const double pi = std::acos(-1.0);
    for (size_t i = 0; i < (order + 1) / 2; ++i) {
        double x = std::cos(pi * (static_cast<double>(i) + 0.75) / (static_cast<double>(order) + 0.5));
        for (size_t iteration = 0; iteration < 100; ++iteration) {
            double p0 = 1.0, p1 = x;
            for (size_t degree = 2; degree <= order; ++degree) {
                double p2 = ((2.0 * degree - 1.0) * x * p1 - (degree - 1.0) * p0) / degree;
                p0 = p1;
                p1 = p2;
            }
            double derivative = static_cast<double>(order) * (p0 - x * p1) / (1.0 - x * x);
            double next = x - p1 / derivative;
            if (std::abs(next - x) < 1e-15) {
                x = next;
                break;
            }
            x = next;
        }
        double p0 = 1.0, p1 = x;
        for (size_t degree = 2; degree <= order; ++degree) {
            double p2 = ((2.0 * degree - 1.0) * x * p1 - (degree - 1.0) * p0) / degree;
            p0 = p1;
            p1 = p2;
        }
        double derivative = static_cast<double>(order) * (p0 - x * p1) / (1.0 - x * x);
        double weight = 2.0 / ((1.0 - x * x) * derivative * derivative);
        nodes[i] = -x;
        nodes[order - 1 - i] = x;
        weights[i] = weight;
        weights[order - 1 - i] = weight;
    }
    return {nodes, weights};
}

std::pair<std::vector<double>, std::vector<double>> hermite(size_t order) {
    if (order == 0)
        throw std::invalid_argument("quadrature order must be positive");
    std::vector<double> nodes(order), weights(order);
    const double pi = std::acos(-1.0);
    const double sqrt_pi = std::sqrt(pi);
    std::vector<std::pair<double, double>> roots;
    roots.reserve(order);
    std::vector<double> positive_roots;
    double x = 0.0;
    for (size_t i = 0; i < (order + 1) / 2; ++i) {
        if (i == 0)
            x = std::sqrt(2.0 * order + 1.0) - 1.85575 * std::pow(2.0 * order + 1.0, -1.0 / 6.0);
        if (i == 1)
            x -= 1.14 * std::pow(static_cast<double>(order), 0.426) / x;
        else if (i == 2)
            x = 1.86 * x - 0.86 * positive_roots[0];
        else if (i == 3)
            x = 1.91 * x - 0.91 * positive_roots[1];
        else if (i > 3)
            x = 2.0 * x - positive_roots[i - 2];
        for (size_t iteration = 0; iteration < 100; ++iteration) {
            double h0 = 1.0, h1 = 2.0 * x;
            for (size_t degree = 2; degree <= order; ++degree) {
                double h2 = 2.0 * x * h1 - 2.0 * (degree - 1.0) * h0;
                h0 = h1;
                h1 = h2;
            }
            double derivative = 2.0 * order * h0;
            double next = x - h1 / derivative;
            if (std::abs(next - x) < 1e-14) {
                x = next;
                break;
            }
            x = next;
        }
        double h0 = 1.0, h1 = 2.0 * x;
        for (size_t degree = 2; degree <= order; ++degree) {
            double h2 = 2.0 * x * h1 - 2.0 * (degree - 1.0) * h0;
            h0 = h1;
            h1 = h2;
        }
        // Store the positive root temporarily in the second slot for the
        // next initial estimate; weights are reconstructed below.
        roots.emplace_back(x, h0);
        positive_roots.push_back(x);
    }
    for (size_t i = 0; i < roots.size(); ++i) {
        const double root = roots[i].first;
        const double h_previous = roots[i].second;
        const double weight = std::pow(2.0, order - 1) * std::tgamma(order + 1) * sqrt_pi /
                              (static_cast<double>(order) * order * h_previous * h_previous);
        nodes[i] = -root;
        weights[i] = weight;
        if (order - 1 - i != i) {
            nodes[order - 1 - i] = root;
            weights[order - 1 - i] = weight;
        }
    }
    return {nodes, weights};
}

void require_positive_finite(double value, const char *name) {
    if (!std::isfinite(value) || value <= 0.0)
        throw py::value_error(std::string(name) + " must be a finite positive scalar");
}

size_t parse_sampling(py::object value) {
    if (py::isinstance<py::bool_>(value) || !py::isinstance<py::int_>(value))
        throw py::value_error("sampling must be a positive integer");
    const auto sampling = value.cast<long long>();
    if (sampling < 1)
        throw py::value_error("sampling must be a positive integer");
    return static_cast<size_t>(sampling);
}

std::vector<double> normalized_weights(py::object values, size_t expected) {
    py::module_ numpy = py::module_::import("numpy");
    if (numpy.attr("iscomplexobj")(values).cast<bool>())
        throw py::value_error("number_weights must be real");
    auto weights = array_like_1d_to_double_vector(std::move(values));
    if (weights.size() != expected || weights.empty())
        throw py::value_error(
            "number_weights must be a nonempty one-dimensional array matching the number of sizes or components");
    double maximum = 0.0;
    for (double weight : weights) {
        if (!std::isfinite(weight) || weight < 0.0)
            throw py::value_error("number_weights must be finite and nonnegative, with at least one positive weight");
        maximum = std::max(maximum, weight);
    }
    if (maximum == 0.0)
        throw py::value_error("number_weights must be finite and nonnegative, with at least one positive weight");
    double total = 0.0;
    for (auto &weight : weights) {
        weight /= maximum;
        total += weight;
    }
    for (auto &weight : weights)
        weight /= total;
    return weights;
}

ParticleSizeDistribution::ParticleSizeDistribution(py::object diameters, py::object number_weights) {
    if (!py::hasattr(diameters, "to"))
        throw py::type_error("diameters must carry length units");
    py::module_ numpy = py::module_::import("numpy");
    if (numpy.attr("iscomplexobj")(diameters.attr("magnitude")).cast<bool>())
        throw py::value_error("diameters must be real");
    diameters_ = quantity_1d_to_meters_vector(std::move(diameters));
    if (diameters_.empty())
        throw py::value_error("diameters must be a nonempty one-dimensional array");
    for (double diameter : diameters_)
        require_positive_finite(diameter, "diameters");
    weights_ = normalized_weights(std::move(number_weights), diameters_.size());
}

ParticleSizeDistribution::ParticleSizeDistribution(std::vector<double> diameters, std::vector<double> weights,
                                                   bool normalize)
    : diameters_(std::move(diameters)), weights_(std::move(weights)) {
    for (double diameter : diameters_)
        require_positive_finite(diameter, "diameters");
    if (normalize) {
        py::gil_scoped_acquire gil;
        weights_ = normalized_weights(py::cast(weights_), diameters_.size());
    }
}

py::object ParticleSizeDistribution::diameters() const {
    return py::cast(diameters_) * get_shared_ureg().attr("meter");
}

py::array_t<double> ParticleSizeDistribution::number_fractions() const {
    py::array_t<double> result(weights_.size());
    std::copy(weights_.begin(), weights_.end(), result.mutable_data());
    return result;
}

ParticleSizeDistribution ParticleSizeDistribution::monodisperse(py::object diameter) {
    double value = quantity_scalar_to_meters(std::move(diameter));
    require_positive_finite(value, "diameter");
    return ParticleSizeDistribution({value}, {1.0}, false);
}

ParticleSizeDistribution ParticleSizeDistribution::uniform(py::object minimum, py::object maximum, size_t sampling) {
    auto [lower, upper] = bounds(std::move(minimum), std::move(maximum));
    auto [nodes, weights] = legendre_unit(sampling);
    for (auto &node : nodes)
        node = lower + (upper - lower) * node;
    return ParticleSizeDistribution(std::move(nodes), std::move(weights), true);
}

ParticleSizeDistribution ParticleSizeDistribution::lognormal(py::object median, double geometric_std, size_t sampling) {
    validate_sampling(sampling);
    py::module_ numpy = py::module_::import("numpy");
    if (!py::hasattr(median, "to") || numpy.attr("asarray")(median.attr("magnitude")).attr("ndim").cast<int>() != 0)
        throw py::value_error("median_diameter must be a finite positive scalar");
    double value = quantity_scalar_to_meters(std::move(median));
    require_positive_finite(value, "median_diameter");
    if (!std::isfinite(geometric_std) || geometric_std < 1.0)
        throw py::value_error("geometric_std must be finite and at least one");
    if (geometric_std == 1.0)
        return ParticleSizeDistribution({value}, {1.0}, false);
    auto [nodes, weights] = hermite(sampling);
    const double factor = std::sqrt(2.0) * std::log(geometric_std);
    for (auto &node : nodes)
        node = value * std::exp(factor * node);
    return ParticleSizeDistribution(std::move(nodes), std::move(weights), true);
}

ParticleSizeDistribution ParticleSizeDistribution::truncated_normal(py::object mean, py::object standard_deviation,
                                                                    py::object minimum, py::object maximum,
                                                                    size_t sampling) {
    double location = quantity_scalar_to_meters(std::move(mean));
    double width = quantity_scalar_to_meters(std::move(standard_deviation));
    require_positive_finite(location, "mean_diameter");
    require_positive_finite(width, "standard_deviation");
    auto [lower, upper] = bounds(std::move(minimum), std::move(maximum));
    auto [nodes, weights] = legendre_unit(sampling);
    double maximum_log_density = -INFINITY;
    for (auto &node : nodes) {
        node = lower + (upper - lower) * node;
        maximum_log_density = std::max(maximum_log_density, -0.5 * std::pow((node - location) / width, 2));
    }
    for (size_t i = 0; i < weights.size(); ++i)
        weights[i] *= std::exp(-0.5 * std::pow((nodes[i] - location) / width, 2) - maximum_log_density);
    return ParticleSizeDistribution(std::move(nodes), std::move(weights), true);
}

ParticleSizeDistribution ParticleSizeDistribution::triangular(py::object minimum, py::object mode, py::object maximum,
                                                              size_t sampling) {
    auto [lower, upper] = bounds(std::move(minimum), std::move(maximum));
    double peak = quantity_scalar_to_meters(std::move(mode));
    require_positive_finite(peak, "mode_diameter");
    if (peak < lower || peak > upper)
        throw py::value_error("mode_diameter must lie between the minimum and maximum diameters");
    auto [unit_nodes, unit_weights] = legendre_unit(sampling);
    std::vector<double> nodes, weights;
    double left_fraction = (peak - lower) / (upper - lower);
    if (peak > lower)
        for (size_t i = 0; i < unit_nodes.size(); ++i) {
            nodes.push_back(lower + (peak - lower) * unit_nodes[i]);
            weights.push_back(2 * unit_weights[i] * unit_nodes[i] * left_fraction);
        }
    if (peak < upper)
        for (size_t i = 0; i < unit_nodes.size(); ++i) {
            nodes.push_back(peak + (upper - peak) * unit_nodes[i]);
            weights.push_back(2 * unit_weights[i] * (1 - unit_nodes[i]) * (1 - left_fraction));
        }
    return ParticleSizeDistribution(std::move(nodes), std::move(weights), true);
}

ParticleSizeDistribution ParticleSizeDistribution::mixture(py::iterable distributions, py::object number_weights) {
    std::vector<ParticleSizeDistribution> components;
    for (py::handle item : distributions) {
        if (!py::isinstance<ParticleSizeDistribution>(item))
            throw py::value_error("distributions must be a nonempty sequence of ParticleSizeDistribution objects");
        components.push_back(item.cast<ParticleSizeDistribution>());
    }
    if (components.empty())
        throw py::value_error("distributions must be a nonempty sequence of ParticleSizeDistribution objects");
    auto component_weights = normalized_weights(std::move(number_weights), components.size());
    std::vector<double> nodes, weights;
    for (size_t i = 0; i < components.size(); ++i)
        for (size_t j = 0; j < components[i].diameters_.size(); ++j) {
            nodes.push_back(components[i].diameters_[j]);
            weights.push_back(component_weights[i] * components[i].weights_[j]);
        }
    return ParticleSizeDistribution(std::move(nodes), std::move(weights), true);
}

std::string ParticleSizeDistribution::repr() const {
    return "ParticleSizeDistribution(samples=" + std::to_string(diameters_.size()) + ", weighting='number')";
}

void ParticleSizeDistribution::validate_sampling(size_t sampling) {
    if (sampling == 0)
        throw py::value_error("sampling must be a positive integer");
}

std::pair<double, double> ParticleSizeDistribution::bounds(py::object minimum, py::object maximum) {
    double lower = quantity_scalar_to_meters(std::move(minimum));
    double upper = quantity_scalar_to_meters(std::move(maximum));
    require_positive_finite(lower, "minimum_diameter");
    require_positive_finite(upper, "maximum_diameter");
    if (lower >= upper)
        throw py::value_error("minimum_diameter must be smaller than maximum_diameter");
    return {lower, upper};
}

std::pair<std::vector<double>, std::vector<double>> ParticleSizeDistribution::legendre_unit(size_t sampling) {
    validate_sampling(sampling);
    auto [nodes, weights] = legendre(sampling);
    for (auto &node : nodes)
        node = (node + 1) / 2;
    for (auto &weight : weights)
        weight /= 2;
    return {nodes, weights};
}
