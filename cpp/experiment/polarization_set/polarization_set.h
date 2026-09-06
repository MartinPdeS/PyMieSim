#pragma once

#include <vector>
#include <complex>
#include <stdexcept>
#include <cmath>
#include <cstddef>

#include <polarization/polarization.h>

using complex128 = std::complex<double>;

class PolarizationSet {
public:
    std::vector<PolarizationState> polarization_states;
    std::vector<double> angles;

    PolarizationSet() = default;

    explicit PolarizationSet(const std::vector<double>& angles_radian)
        : angles(angles_radian)
    {
        if (angles_radian.empty()) {
            throw std::invalid_argument("Polarization angles cannot be empty.");
        }

        this->polarization_states.reserve(angles_radian.size());

        for (const double angle : angles_radian) {
            if (!std::isfinite(angle)) {
                throw std::invalid_argument("Polarization angles must be finite values.");
            }

            this->polarization_states.emplace_back(angle);
        }
    }

    explicit PolarizationSet(const PolarizationState& polarization_state)
    {
        this->polarization_states.push_back(polarization_state);

        if (!polarization_state.is_empty) {
            if (polarization_state.jones_vector.size() != 2) {
                throw std::invalid_argument(
                    "PolarizationState must have exactly 2 Jones vector components [Ex, Ey], unless it is empty."
                );
            }

            if (std::isfinite(polarization_state.angle)) {
                this->angles.push_back(polarization_state.angle);
            }
        }
    }

    explicit PolarizationSet(const double angle_radian)
        : angles({angle_radian})
    {
        if (!std::isfinite(angle_radian)) {
            throw std::invalid_argument("Polarization angle must be a finite value.");
        }

        this->polarization_states.emplace_back(angle_radian);
    }

    explicit PolarizationSet(const std::vector<std::vector<complex128>>& jones_vectors)
    {
        if (jones_vectors.empty()) {
            throw std::invalid_argument("Jones vectors cannot be empty.");
        }

        this->polarization_states.reserve(jones_vectors.size());

        for (const std::vector<complex128>& jones_vector : jones_vectors) {
            if (jones_vector.size() != 2) {
                throw std::invalid_argument(
                    "Each Jones vector must have exactly 2 components [Ex, Ey]."
                );
            }

            this->polarization_states.emplace_back(jones_vector);
        }
    }

    explicit PolarizationSet(const std::vector<PolarizationState>& polarization_states)
        : polarization_states(polarization_states)
    {
        if (polarization_states.empty()) {
            throw std::invalid_argument("Polarization states cannot be empty.");
        }

        this->angles.clear();
        this->angles.reserve(polarization_states.size());

        for (const PolarizationState& state : polarization_states) {
            if (state.is_empty) {
                continue;
            }

            if (state.jones_vector.size() != 2) {
                throw std::invalid_argument(
                    "Each PolarizationState must have exactly 2 Jones vector components [Ex, Ey], unless it is empty."
                );
            }

            if (std::isfinite(state.angle)) {
                this->angles.push_back(state.angle);
            }
        }
    }

    size_t size() const noexcept
    {
        return this->polarization_states.size();
    }

    size_t number_of_states() const noexcept
    {
        return this->polarization_states.size();
    }

    PolarizationState operator[](size_t index) const
    {
        if (index >= this->polarization_states.size()) {
            throw std::out_of_range("Polarization index out of range.");
        }

        return this->polarization_states[index];
    }
};