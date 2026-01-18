#pragma once

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <fibonacci/fibonacci.h>

/**
 * @brief Snell and Fresnel utilities for a planar interface between two real index media.
 *
 * This class centralizes three distinct effects that you may want to enable or disable:
 *
 * 1) Angle mapping (Snell):
 *    Maps detector NA limits (defined in detector medium) into angular bounds in the scatterer medium.
 *
 * 2) Angle clamp:
 *    Enforces physical bounds on angles (sin <= 1) and optionally the TIR critical angle clamp when n_s > n_d.
 *    If you disable clamping, the class will NOT clamp to pi/2 or to the critical angle. In that case,
 *    you must decide how to handle invalid asin arguments and TIR in the caller. Here we choose to throw.
 *
 * 3) Fresnel transmission:
 *    Computes Fresnel amplitude transmission coefficients t_s and t_p on a sampling mesh.
 *    If disabled, outputs are filled with 1.0 for forward directions and 0.0 for backward directions.
 *
 * Notes:
 * - Indices are assumed real and strictly positive.
 * - Fresnel coefficients are defined for plane-wave incidence from scatterer medium (n_s) to detector medium (n_d).
 */
class SnellFresnelInterface {
public:
    /**
     * @brief Construct an interface with default indices n_s = 1 and n_d = 1.
     *
     * The indices should typically be set explicitly using set_media or the two argument constructor
     * before calling any computation method.
     */
    SnellFresnelInterface() = default;

    /**
     * @brief Construct an interface with explicit propagation indices.
     *
     * @param scatterer_medium_refractive_index Refractive index n_s of the incident medium.
     * @param detector_medium_refractive_index  Refractive index n_d of the transmitted medium.
     *
     * @throws std::invalid_argument if any refractive index is not strictly positive.
     */
    SnellFresnelInterface(
        double scatterer_medium_refractive_index,
        double detector_medium_refractive_index)
    {
        this->set_media(scatterer_medium_refractive_index, detector_medium_refractive_index);
    }

    /**
     * @brief Set the propagation indices of the two media defining the interface.
     *
     * @param scatterer_medium_refractive_index Refractive index n_s of the incident medium.
     * @param detector_medium_refractive_index  Refractive index n_d of the transmitted medium.
     *
     * @throws std::invalid_argument if any refractive index is not strictly positive.
     */
    void set_media(
        double scatterer_medium_refractive_index,
        double detector_medium_refractive_index)
    {
        if (scatterer_medium_refractive_index <= 0.0 || detector_medium_refractive_index <= 0.0) {
            throw std::invalid_argument("refractive indices must be positive.");
        }

        this->scatterer_medium_refractive_index = scatterer_medium_refractive_index;
        this->detector_medium_refractive_index = detector_medium_refractive_index;
    }

    /**
     * @brief Enable or disable Snell angle mapping for NA -> scatterer angle bounds.
     *
     * If disabled, compute_snell_mapped_angle_bounds will treat the detector medium as the same as the
     * scatterer medium for bounds purposes (i.e. theta_scatterer = theta_detector).
     */
    void set_enable_angle_mapping(bool enabled) { this->enable_angle_mapping = enabled; }

    /**
     * @brief Enable or disable angle clamping (pi/2 and critical-angle clamp).
     *
     * If disabled, the class will NOT clamp asin arguments and will throw if a physically invalid value
     * is encountered (e.g., NA/nd > 1 or Snell produces sin(theta) > 1).
     */
    void set_enable_angle_clamp(bool enabled) { this->enable_angle_clamp = enabled; }

    /**
     * @brief Enable or disable Fresnel transmission coefficients computation.
     *
     * If disabled, compute_fresnel_amplitude_transmission_on_mesh returns:
     * - t_s = 1, t_p = 1 for forward directions (cos(theta_s) > 0)
     * - t_s = 0, t_p = 0 for backward directions (cos(theta_s) <= 0)
     *
     * This keeps the "forward hemisphere only" behavior without applying any interface loss.
     */
    void set_enable_fresnel(bool enabled) { this->enable_fresnel = enabled; }

    /**
     * @brief Get current toggle states.
     */
    bool get_enable_angle_mapping() const { return this->enable_angle_mapping; }
    bool get_enable_angle_clamp() const { return this->enable_angle_clamp; }
    bool get_enable_fresnel() const { return this->enable_fresnel; }

    /**
     * @brief Compute angular bounds in the scatterer medium that correspond to detector NA limits.
     *
     * The detector numerical aperture is interpreted in the detector medium:
     *     NA_d = n_d * sin(theta_d)
     * so:
     *     theta_d = asin( NA_d / n_d ).
     *
     * If angle mapping is enabled, the corresponding incident angle in the scatterer medium is obtained by Snell law:
     *     n_s * sin(theta_s) = n_d * sin(theta_d)
     * so:
     *     sin(theta_s) = (n_d / n_s) * sin(theta_d).
     *
     * If angle mapping is disabled, we take:
     *     theta_s = theta_d
     * i.e. we ignore refraction and interpret the NA in the same medium.
     *
     * If angle clamp is enabled:
     * - asin arguments are clamped into [-1, 1]
     * - if n_s > n_d and angle mapping is enabled, theta_s is clamped to the critical angle:
     *     theta_c = asin(n_d / n_s)
     *
     * If angle clamp is disabled:
     * - this method throws if any asin argument is out of range
     * - this method does NOT apply critical angle clamping
     *
     * @param detector_numerical_aperture        Detector NA upper bound (dimensionless).
     * @param detector_cache_numerical_aperture  Detector cache NA lower bound (dimensionless).
     * @param[out] theta_scatterer_max           Maximum accepted polar angle in scatterer medium (radians).
     * @param[out] theta_scatterer_min           Minimum accepted polar angle in scatterer medium (radians).
     *
     * @throws std::invalid_argument if the mapped min bound exceeds the mapped max bound.
     * @throws std::domain_error if angle clamp is disabled and invalid values occur.
     */
    void compute_snell_mapped_angle_bounds(
        double detector_numerical_aperture,
        double detector_cache_numerical_aperture,
        double& theta_scatterer_max,
        double& theta_scatterer_min) const
    {
        const double ns = this->scatterer_medium_refractive_index;
        const double nd = this->detector_medium_refractive_index;

        const double theta_detector_max = this->detector_na_to_theta_in_detector_medium(detector_numerical_aperture);
        const double theta_detector_min = this->detector_na_to_theta_in_detector_medium(detector_cache_numerical_aperture);

        if (!this->enable_angle_mapping) {
            theta_scatterer_max = theta_detector_max;
            theta_scatterer_min = theta_detector_min;
        } else {
            const double sin_theta_scatterer_max = (nd / ns) * std::sin(theta_detector_max);
            const double sin_theta_scatterer_min = (nd / ns) * std::sin(theta_detector_min);

            theta_scatterer_max = this->asin_unit(sin_theta_scatterer_max);
            theta_scatterer_min = this->asin_unit(sin_theta_scatterer_min);

            if (this->enable_angle_clamp && (ns > nd)) {
                const double theta_critical = this->asin_unit(nd / ns);
                theta_scatterer_max = std::min(theta_scatterer_max, theta_critical);
                theta_scatterer_min = std::min(theta_scatterer_min, theta_critical);
            }
        }

        if (theta_scatterer_max < theta_scatterer_min) {
            throw std::invalid_argument("Mapped cache NA cannot be larger than mapped detector NA.");
        }
    }

    /**
     * @brief Compute Fresnel amplitude transmission coefficients on a FibonacciMesh.
     *
     * If Fresnel is enabled:
     * - Computes t_s and t_p for incidence from n_s to n_d.
     * - Uses Snell to determine theta_d; if sin(theta_d) >= 1 then TIR -> t_s = t_p = 0.
     *
     * If Fresnel is disabled:
     * - Sets t_s = t_p = 1 for forward directions (cos(theta_s) > 0),
     * - Sets t_s = t_p = 0 for backward directions (cos(theta_s) <= 0).
     *
     * Angle clamp interaction:
     * - This method always treats sin(theta_d) >= 1 as TIR and sets coefficients to 0, because
     *   Fresnel transmission is physically undefined there for propagating transmission.
     * - Disabling "angle clamp" does not disable this TIR check, because it is part of Fresnel physics.
     *
     * @param fibonacci_mesh The direction sampling mesh used by the detector.
     * @param[out] interface_t_s Output vector of t_s values per sample (size = mesh sampling).
     * @param[out] interface_t_p Output vector of t_p values per sample (size = mesh sampling).
     */
    void compute_fresnel_amplitude_transmission_on_mesh(
        const FibonacciMesh& fibonacci_mesh,
        std::vector<double>& interface_t_s,
        std::vector<double>& interface_t_p) const
    {
        const double ns = this->scatterer_medium_refractive_index;
        const double nd = this->detector_medium_refractive_index;

        const size_t sampling = fibonacci_mesh.cartesian.z.size();

        interface_t_s.assign(sampling, 0.0);
        interface_t_p.assign(sampling, 0.0);

        for (size_t i = 0; i < sampling; ++i) {

            const double cos_theta_s = clamp_m1_p1(fibonacci_mesh.cartesian.z[i]);

            if (cos_theta_s <= 0.0) {
                interface_t_s[i] = 0.0;
                interface_t_p[i] = 0.0;
                continue;
            }

            if (!this->enable_fresnel) {
                interface_t_s[i] = 1.0;
                interface_t_p[i] = 1.0;
                continue;
            }

            const double sin_theta_s_sq = std::max(0.0, 1.0 - cos_theta_s * cos_theta_s);
            const double sin_theta_s = std::sqrt(sin_theta_s_sq);

            // Snell: n_s sin(theta_s) = n_d sin(theta_d) -> sin(theta_d) = (n_s/n_d) sin(theta_s)
            const double sin_theta_d = (ns / nd) * sin_theta_s;

            // TIR: evanescent transmitted wave => no amplitude transmission for propagating field
            if (sin_theta_d >= 1.0) {
                interface_t_s[i] = 0.0;
                interface_t_p[i] = 0.0;
                continue;
            }

            const double cos_theta_d = std::sqrt(std::max(0.0, 1.0 - sin_theta_d * sin_theta_d));

            // Fresnel amplitude transmission coefficients (E field)
            const double denom_s = ns * cos_theta_s + nd * cos_theta_d;
            const double denom_p = nd * cos_theta_s + ns * cos_theta_d;

            if (denom_s == 0.0 || denom_p == 0.0) {
                interface_t_s[i] = 0.0;
                interface_t_p[i] = 0.0;
                continue;
            }

            interface_t_s[i] = (2.0 * ns * cos_theta_s) / denom_s;
            interface_t_p[i] = (2.0 * ns * cos_theta_s) / denom_p;
        }
    }

private:
    /**
     * @brief Refractive index of the incident propagation medium n_s.
     */
    double scatterer_medium_refractive_index = 1.0;

    /**
     * @brief Refractive index of the transmitted propagation medium n_d.
     */
    double detector_medium_refractive_index = 1.0;

    /**
     * @brief Toggle Snell mapping of angles.
     */
    bool enable_angle_mapping = true;

    /**
     * @brief Toggle clamping for asin arguments and critical-angle clamping.
     */
    bool enable_angle_clamp = true;

    /**
     * @brief Toggle Fresnel transmission coefficients.
     */
    bool enable_fresnel = true;

private:
    /**
     * @brief Clamp a value into [-1, 1] for robust trigonometric operations.
     */
    static inline double clamp_m1_p1(double x) {
        return std::max(-1.0, std::min(1.0, x));
    }

    /**
     * @brief asin(x) with configurable clamping behavior.
     *
     * If enable_angle_clamp is true, x is clamped to [-1, 1] before asin.
     * If enable_angle_clamp is false, throws if x is outside [-1, 1].
     */
    double asin_unit(double x) const
    {
        if (this->enable_angle_clamp) {
            return std::asin(clamp_m1_p1(x));
        }

        if (x < -1.0 || x > 1.0) {
            throw std::domain_error("asin argument out of range and angle clamp disabled.");
        }
        return std::asin(x);
    }

    /**
     * @brief Convert a detector NA into the corresponding acceptance half angle in the detector medium.
     *
     * Interpretation:
     *     NA_d = n_d * sin(theta_d)  ->  theta_d = asin( NA_d / n_d ).
     *
     * If enable_angle_clamp is true, NA_d / n_d is clamped to [-1, 1].
     * If enable_angle_clamp is false, this throws when NA_d / n_d is out of [-1, 1].
     *
     * @param detector_numerical_aperture Detector NA (dimensionless).
     * @return Acceptance half angle theta_d in radians, within [0, pi/2] if clamping is enabled.
     */
    double detector_na_to_theta_in_detector_medium(double detector_numerical_aperture) const
    {
        const double nd = this->detector_medium_refractive_index;
        const double arg = detector_numerical_aperture / nd;

        return this->asin_unit(arg);
    }
};
