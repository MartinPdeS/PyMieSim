#pragma once

#include <single/detector/base.h>


class CoherentMode : public BaseDetector {
public:
    std::string mode_number;
    double rotation = 0.0;

public:
    CoherentMode() = default;

    CoherentMode(
        const std::string &_mode_number,
        const size_t _sampling,
        const double _numerical_aperture,
        const double _cache_numerical_aperture,
        const double _phi_offset,
        const double _gamma_offset,
        const PolarizationState& _polarization_filter,
        const double _rotation,
        const bool _mean_coupling,
        const std::shared_ptr<BaseMedium> _medium
    )
    :   BaseDetector(
            _sampling,
            _numerical_aperture,
            _cache_numerical_aperture,
            _phi_offset,
            _gamma_offset,
            _polarization_filter,
            std::move(_medium),
            _mean_coupling
        ),
        mode_number(_mode_number),
        rotation(_rotation)
    {
        this->parse_mode(mode_number);
        this->mode_field = ModeField(this->mode_id);
    }

    CoherentMode(
        const std::string &_mode_number,
        const size_t _sampling,
        const double _numerical_aperture,
        const double _cache_numerical_aperture,
        const double _phi_offset,
        const double _gamma_offset,
        const PolarizationState& _polarization_filter,
        const double _rotation,
        const bool _mean_coupling,
        const double _medium
    )
    :   CoherentMode(
            _mode_number,
            _sampling,
            _numerical_aperture,
            _cache_numerical_aperture,
            _phi_offset,
            _gamma_offset,
            _polarization_filter,
            _rotation,
            _mean_coupling,
            std::make_shared<ConstantMedium>(_medium)
        )
    {}

    double get_coupling(std::shared_ptr<BaseScatterer> scatterer, std::shared_ptr<BaseSource> source) override;

    [[nodiscard]] std::vector<complex128> get_structured_scalarfield(const size_t sampling) const override;

    /**
     * @brief Initializes the detector's mesh based on the scatterer's properties and the interface conditions.
     * @param scatterer The scatterer for which the mesh is initialized.
     * @param rotation The rotation angle to apply to the mesh.
     * @note This function sets up the Fibonacci mesh for the detector, computes the scalar field based on the mode field, and applies the Fresnel transmission coefficients to the fields.
     */
    void initialize_mesh(const std::shared_ptr<BaseScatterer> scatterer) override;

private:

    /**
     * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    double get_coupling_point(std::shared_ptr<BaseScatterer> scatterer, std::shared_ptr<BaseSource> source);

    /**
     * @brief Computes the coupling coefficient for a given scatterer, considering coherent modes.
     * @param scatterer The scatterer for which the coupling coefficient is computed.
     * @return The coupling coefficient.
     * @note This function computes the coupling coefficient based on the mode field and the scatterer's properties.
     */
    double get_coupling_mean(std::shared_ptr<BaseScatterer> scatterer, std::shared_ptr<BaseSource> source);

    /**
     * @brief Computes the projected fields based on the theta and phi fields.
     * @param theta_field The field in the theta direction.
     * @param phi_field The field in the phi direction.
     * @return A tuple containing the horizontal and vertical projections of the fields.
     * @note This function computes the projections based on the spherical coordinates of the Fibonacci mesh.
     */
    std::tuple<std::vector<complex128>, std::vector<complex128>>
    get_projected_farfields(const std::vector<complex128>& theta_field, const std::vector<complex128>& phi_field) const;

    /**
     * @brief Applies a scalar field to the horizontal and vertical projections of the coupling coefficients.
     * @param field0 The horizontal projection field.
     * @param field1 The vertical projection field.
     * @note This function modifies the input fields in place, applying the scalar field to both projections.
     */
    void apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const;

    /**
     * @brief Applies a polarization filter to the coupling coefficients.
     * @tparam T Type of the coupling coefficients.
     * @param coupling_theta The coupling coefficient for the theta direction.
     * @param coupling_phi The coupling coefficient for the phi direction.
     * @param polarization_filter The polarization filter value.
     */
    template <typename T> inline void apply_polarization_filter(T& coupling_theta, T& coupling_phi) const;
};

