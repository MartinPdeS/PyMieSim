#pragma once

#include <array>
#include <complex>
#include <cstddef>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <limits>


using complex128 = std::complex<double>;

class PolarizationState {
public:
    std::vector<complex128> jones_vector;
    double angle = std::numeric_limits<double>::quiet_NaN();

    bool is_empty = false;

    PolarizationState() {
        this->is_empty = true;
    };

    PolarizationState(const std::vector<complex128>& jones_vector)
        : jones_vector(jones_vector)
    {
        if (this->jones_vector.size() != 2) {
            throw std::invalid_argument("jones_vector must have exactly 2 components [Ex, Ey].");
        }
        this->is_empty = false;
    }

    PolarizationState(const double angle_radian)
    {
        this->jones_vector = {std::cos(angle_radian), std::sin(angle_radian)};
        this->angle = angle_radian;
        this->is_empty = false;
    }
};


class LeftCircular : public PolarizationState {
public:
    LeftCircular() : PolarizationState({1.0 / std::sqrt(2), 1.0 / std::sqrt(2) * complex128(0, 1)}) {}
};


class RightCircular : public PolarizationState {
public:
    RightCircular() : PolarizationState({1.0 / std::sqrt(2), -1.0 / std::sqrt(2) * complex128(0, 1)}) {}
};