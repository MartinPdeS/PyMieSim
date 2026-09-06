#pragma once

#include <single/detector/base.h>

class IntegratingSphere : public BaseDetector {
public:
    /**
     * @brief Construct an integrating sphere detector that collects over 4π steradians.
     *
     * This detector samples directions over the full sphere (theta in [0, pi]).
     * It is intended as a power collector: it integrates the scattered field intensity
     * over all directions, optionally applying a transmission model at a planar interface
     * between scatterer medium and detector medium (SnellFresnelInterface).
     *
     * Notes:
     * - numerical_aperture and cache_numerical_aperture are ignored for this detector.
     * - phi_offset and gamma_offset can still be used to rotate the sampling pattern,
     *   but they do not restrict the acceptance.
     * - If you want an ideal integrating sphere with no interface physics, disable
     *   Fresnel transmission in SnellFresnelInterface (or set n_s = n_d and disable).
     */

    IntegratingSphere(
        const size_t sampling,
        const PolarizationState& polarization_filter)
    :   BaseDetector(
            sampling,
            0.0,  /* numerical_aperture */
            0.0,  /* cache_numerical_aperture */
            0.0,  /* phi_offset */
            0.0,  /* gamma_offset */
            polarization_filter,
            std::make_shared<ConstantMedium>(1.0),  /* medium */
            false
        )
    {
        this->medium = std::make_shared<ConstantMedium>(1.0);
        this->parse_mode("NC00");
        this->mode_field = ModeField(this->mode_id);
    }

    /**
     * @brief Compute the collected power from the scatterer over 4π.
     *
     * For incoherent detection, the coupling is proportional to the integral
     * of |E_theta|^2 + |E_phi|^2 over the full solid angle.
     *
     * If Fresnel transmission is enabled in snell_interface, a scalar amplitude
     * transmission is applied to the fields before intensity computation.
     */
    double get_coupling(std::shared_ptr<BaseScatterer> scatterer, std::shared_ptr<BaseSource> source) override;

    /**
     * @brief Not used for integrating sphere.
     *
     * This detector does not have a mode field. The function returns an empty vector.
     */
    [[nodiscard]] std::vector<complex128> get_structured_scalarfield(const size_t sampling) const override;


    void initialize_mesh(const std::shared_ptr<BaseScatterer> scatterer) override;

};
