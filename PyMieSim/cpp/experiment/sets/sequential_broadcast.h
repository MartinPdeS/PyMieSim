#pragma once

#include <pybind11/pybind11.h>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <sstream>

namespace py = pybind11;

inline size_t resolve_target_size(
    const std::optional<size_t>& total_size,
    const std::vector<std::pair<std::string, std::vector<double>>>& named_vectors
) {
    size_t max_size = 0;
    for (const auto& kv : named_vectors) {
        max_size = std::max(max_size, kv.second.size());
    }

    if (total_size.has_value()) {
        if (*total_size == 0) {
            throw std::invalid_argument("total_size must be a positive integer.");
        }
        if (*total_size < max_size) {
            std::ostringstream oss;
            oss << "total_size (" << *total_size << ") is smaller than the maximum input size (" << max_size << ").";
            throw std::invalid_argument(oss.str());
        }
        return *total_size;
    }

    return std::max<size_t>(1, max_size);
}

inline std::vector<double> broadcast_vector_double(
    const std::string& name,
    const std::vector<double>& values,
    const size_t target_size
) {
    if (values.empty()) {
        std::ostringstream oss;
        oss << "Parameter '" << name << "' is empty. Provide a scalar or a non empty array.";
        throw std::invalid_argument(oss.str());
    }

    if (values.size() == 1) {
        return std::vector<double>(target_size, values[0]);
    }

    if (values.size() != target_size) {
        std::ostringstream oss;
        oss << "Inconsistent sizes: '" << name << "' has size " << values.size()
            << " but expected 1 or " << target_size << ".";
        throw std::invalid_argument(oss.str());
    }

    return values;
}

inline std::vector<unsigned> broadcast_vector_unsigned(
    const std::string& name,
    const std::vector<unsigned>& values,
    const size_t target_size
) {
    if (values.empty()) {
        std::ostringstream oss;
        oss << "Parameter '" << name << "' is empty. Provide a scalar or a non empty array.";
        throw std::invalid_argument(oss.str());
    }

    if (values.size() == 1) {
        return std::vector<unsigned>(target_size, values[0]);
    }

    if (values.size() != target_size) {
        std::ostringstream oss;
        oss << "Inconsistent sizes: '" << name << "' has size " << values.size()
            << " but expected 1 or " << target_size << ".";
        throw std::invalid_argument(oss.str());
    }

    return values;
}

inline std::vector<std::string> broadcast_vector_string(
    const std::string& name,
    const std::vector<std::string>& values,
    const size_t target_size
) {
    if (values.empty()) {
        std::ostringstream oss;
        oss << "Parameter '" << name << "' is empty. Provide a scalar or a non empty array.";
        throw std::invalid_argument(oss.str());
    }

    if (values.size() == 1) {
        return std::vector<std::string>(target_size, values[0]);
    }

    if (values.size() != target_size) {
        std::ostringstream oss;
        oss << "Inconsistent sizes: '" << name << "' has size " << values.size()
            << " but expected 1 or " << target_size << ".";
        throw std::invalid_argument(oss.str());
    }

    return values;
}

inline std::optional<size_t> parse_optional_total_size(const py::object& total_size_obj) {
    if (total_size_obj.is_none()) {
        return std::nullopt;
    }

    if (!py::isinstance<py::int_>(total_size_obj)) {
        throw std::invalid_argument("total_size must be an int or None.");
    }

    const auto value = total_size_obj.cast<long long>();
    if (value <= 0) {
        throw std::invalid_argument("total_size must be a positive integer.");
    }

    return static_cast<size_t>(value);
}

inline size_t resolve_target_size_from_sizes(
    const std::optional<size_t>& total_size,
    const std::vector<std::pair<std::string, size_t>>& named_sizes
) {
    size_t max_size = 0;
    for (const auto& kv : named_sizes) {
        max_size = std::max(max_size, kv.second);
    }

    if (total_size.has_value()) {
        if (*total_size == 0) {
            throw std::invalid_argument("total_size must be a positive integer.");
        }
        if (*total_size < max_size) {
            std::ostringstream oss;
            oss << "total_size (" << *total_size << ") is smaller than the maximum input size (" << max_size << ").";
            throw std::invalid_argument(oss.str());
        }
        return *total_size;
    }

    return std::max<size_t>(1, max_size);
}