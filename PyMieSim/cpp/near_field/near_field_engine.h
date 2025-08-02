#pragma once

#include <complex>
#include <vector>
#include <tuple>
#include <cmath>
#include "scatterer/base_scatterer/base_scatterer.h"
#include "utils/special_function.cpp"

using complex128 = std::complex<double>;

/**
 * @brief Near-field electromagnetic computation engine for spherical scatterers.
 *
 * This class implements near-field computation using Mie scattering theory.
 * The electromagnetic fields are computed using the multipole expansion
 * with vector spherical harmonics (VSH).
 *
 * For external fields (r > a):
 * E(r) = Σ [a_n M_n^(3) + b_n N_n^(3)]
 *
 * For internal fields (r < a):
 * E(r) = Σ [c_n M_n^(1) + d_n N_n^(1)]
 *
 * Where M_n and N_n are vector spherical harmonics, and superscripts
 * indicate the type of spherical Bessel function used.
 *
 * @note Currently supports only spherical geometries. Cylinder support
 *       requires implementation of cn/dn coefficients for infinite cylinders.
 */
class NearFieldEngine
{
public:
    /**
     * @brief Construct near-field computation engine.
     *
     * @param scatterer Reference to the scatterer object containing
     *                  geometry and material properties.
     * @param max_order Maximum multipole order for field computation.
     *                  If 0, automatically determined from size parameter.
     */
    explicit NearFieldEngine(const BaseScatterer& scatterer, size_t max_order = 0);

    /**
     * @brief Compute electromagnetic fields at specified points.
     *
     * @param points 3D coordinates where fields are computed, shape (N, 3).
     * @param field_components List of field components to compute:
     *                        ["Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "|E|", "|H|", "|E|²", "|H|²"]
     * @param normalize_fields Whether to normalize by incident field amplitude.
     *
     * @return Computed field values organized by component and spatial point.
     *
     * @throws std::invalid_argument If points are invalid or field components unknown.
     * @throws std::runtime_error If scatterer doesn't support near-field computation.
     */
    std::map<std::string, std::vector<complex128>> compute_fields(
        const std::vector<std::array<double, 3>>& points,
        const std::vector<std::string>& field_components,
        bool normalize_fields = true
    ) const;

    /**
     * @brief Compute electric field components at a single point.
     *
     * @param x,y,z Cartesian coordinates of observation point.
     * @return Tuple of (Ex, Ey, Ez) complex field components.
     */
    std::tuple<complex128, complex128, complex128> compute_electric_field(
        double x, double y, double z
    ) const;

    /**
     * @brief Compute magnetic field components at a single point.
     *
     * @param x,y,z Cartesian coordinates of observation point.
     * @return Tuple of (Hx, Hy, Hz) complex field components.
     */
    std::tuple<complex128, complex128, complex128> compute_magnetic_field(
        double x, double y, double z
    ) const;

    /**
     * @brief Check if computation point is inside the scatterer.
     *
     * @param x,y,z Cartesian coordinates.
     * @return True if point is inside scatterer boundary.
     */
    bool is_point_inside_scatterer(double x, double y, double z) const;

    /**
     * @brief Get field enhancement factor at a point.
     *
     * Field enhancement is defined as |E_total|/|E_incident| for electric field
     * or |H_total|/|H_incident| for magnetic field.
     *
     * @param x,y,z Cartesian coordinates.
     * @param field_type "electric" or "magnetic".
     * @return Field enhancement factor (real-valued).
     */
    double get_field_enhancement(double x, double y, double z, const std::string& field_type) const;

private:
    const BaseScatterer& m_scatterer;  ///< Reference to scatterer object
    size_t m_max_order;                ///< Maximum multipole order
    double m_radius;                   ///< Scatterer radius
    complex128 m_refractive_index;     ///< Scatterer refractive index
    double m_medium_index;             ///< Medium refractive index
    double m_wavenumber;               ///< Wavenumber in medium

    /**
     * @brief Convert Cartesian to spherical coordinates.
     *
     * @param x,y,z Cartesian coordinates.
     * @return Tuple of (r, theta, phi) in spherical coordinates.
     */
    std::tuple<double, double, double> cartesian_to_spherical(double x, double y, double z) const;

    /**
     * @brief Compute vector spherical harmonics M_mn and N_mn.
     *
     * These are the fundamental building blocks of electromagnetic multipole fields.
     *
     * @param n Multipole order (n >= 1).
     * @param m Azimuthal index (-n <= m <= n).
     * @param r,theta,phi Spherical coordinates.
     * @param bessel_type 1 for j_n (regular), 3 for h_n^(1) (outgoing).
     *
     * @return Tuple of vector components for M_mn and N_mn in spherical coordinates.
     */
    std::tuple<
        std::tuple<complex128, complex128, complex128>,  // M_mn = (M_r, M_theta, M_phi)
        std::tuple<complex128, complex128, complex128>   // N_mn = (N_r, N_theta, N_phi)
    > compute_vector_spherical_harmonics(
        int n, int m, double r, double theta, double phi, int bessel_type
    ) const;

    /**
     * @brief Transform vector components from spherical to Cartesian coordinates.
     *
     * @param v_r,v_theta,v_phi Vector components in spherical coordinates.
     * @param theta,phi Angular coordinates for transformation.
     *
     * @return Tuple of (v_x, v_y, v_z) in Cartesian coordinates.
     */
    std::tuple<complex128, complex128, complex128> spherical_to_cartesian_vector(
        complex128 v_r, complex128 v_theta, complex128 v_phi,
        double theta, double phi
    ) const;

    /**
     * @brief Compute incident field at observation point.
     *
     * This provides the reference field for normalization and enhancement calculations.
     *
     * @param x,y,z Cartesian coordinates.
     * @param field_type "electric" or "magnetic".
     *
     * @return Incident field vector (Ex, Ey, Ez) or (Hx, Hy, Hz).
     */
    std::tuple<complex128, complex128, complex128> compute_incident_field(
        double x, double y, double z, const std::string& field_type
    ) const;

    /**
     * @brief Validate that scatterer supports near-field computation.
     *
     * @throws std::runtime_error If scatterer type is not supported.
     */
    void validate_scatterer_compatibility() const;
};
