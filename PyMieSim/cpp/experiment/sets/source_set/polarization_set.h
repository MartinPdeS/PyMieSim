#pragma once

#include <single/polarization/polarization.h>

using complex128 = std::complex<double>;



class PolarizationSet {
public:
    std::vector<PolarizationState> polarization_states;
    std::vector<double> angles;

    PolarizationSet() = default;

    explicit PolarizationSet(const std::vector<double>& angles_radian)
        : angles(angles_radian)
    {
        for (double angle : angles) {
            // jones_vector.push_back({std::cos(angle), std::sin(angle)});
            PolarizationState temp(angle);
            polarization_states.emplace_back(temp);
            if (!std::isfinite(angle)) {
                throw std::invalid_argument("Polarization angles must be finite values.");
            }
        }
        if (angles_radian.empty()) {
            throw std::invalid_argument("Polarization angles cannot be empty.");
        }
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

