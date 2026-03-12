#pragma once

#include <polarization/polarization.h>

using complex128 = std::complex<double>;



class PolarizationSet {
public:
    std::vector<PolarizationState> polarization_states;
    std::vector<double> angles;

    PolarizationSet() = default;

    /**
     * Initializes the PolarizationSet with a vector of angles in radians. Each angle corresponds to a polarization state, where 0 radians corresponds to linear polarization along the x-axis, and π/4 radians corresponds to circular polarization.
     * @param angles_radian A vector of angles in radians representing the orientation of the polarization states. Each angle must be a finite value, and the vector cannot be empty.
     * @throws std::invalid_argument if the input vector is empty or contains non-finite values.
     */
    explicit PolarizationSet(const std::vector<double>& angles_radian)
        : angles(angles_radian)
    {
        if (angles_radian.empty()) {
            throw std::invalid_argument("Polarization angles cannot be empty.");
        }
        for (double angle : angles) {
            if (!std::isfinite(angle)) {
                throw std::invalid_argument("Polarization angles must be finite values.");
            }
            polarization_states.emplace_back(angle);

        }

    }

    /**
     * Initializes the PolarizationSet with a single angle in radians. This constructor creates a PolarizationSet with one polarization state corresponding to the provided angle.
     * @param angle_radian An angle in radians representing the orientation of the polarization state. The angle must be a finite value.
     * @throws std::invalid_argument if the input angle is not a finite value.
     */
    explicit PolarizationSet(const double angle_radian)
        : angles({angle_radian})
    {
        if (!std::isfinite(angle_radian)) {
            throw std::invalid_argument("Polarization angle must be a finite value.");
        }
        polarization_states.emplace_back(angle_radian);
    }

    /**
     * Initializes the PolarizationSet with a vector of Jones vectors. Each Jones vector must have exactly 2 components [Ex, Ey], representing the x and y components of the electric field, respectively.
     * @param jones_vectors A vector of Jones vectors, where each Jones vector is a vector of two complex numbers [Ex, Ey].
     */
    explicit PolarizationSet(const std::vector<std::vector<complex128>>& jones_vectors)
    {
        for (const std::vector<complex128>& jones_vector : jones_vectors) {
            if (jones_vector.size() != 2) {
                throw std::invalid_argument("Each Jones vector must have exactly 2 components [Ex, Ey].");
            }
            PolarizationState temp(jones_vector);
            polarization_states.emplace_back(temp);
        }
    }


    explicit PolarizationSet(const std::vector<PolarizationState>& polarization_states)
        : polarization_states(polarization_states)
    {
        if (polarization_states.empty()) {
            throw std::invalid_argument("Polarization states cannot be empty.");
        }
        for (const PolarizationState& state : polarization_states) {
            if (state.jones_vector.size() != 2) {
                throw std::invalid_argument("Each PolarizationState must have a Jones vector with exactly 2 components [Ex, Ey].");
            }
        }
    }

    /**
     * Returns the number of polarization states in the set.
     * @return The number of polarization states.
     */
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

