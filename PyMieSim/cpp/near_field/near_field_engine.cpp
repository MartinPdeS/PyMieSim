#include "near_field_engine.h"
#include <stdexcept>
#include <map>
#include <cmath>

NearFieldEngine::NearFieldEngine(const BaseScatterer& scatterer, size_t max_order)
    : m_scatterer(scatterer)
    , m_max_order(max_order == 0 ? scatterer.max_order : max_order)
    , m_radius(scatterer.size_parameter / scatterer.medium_refractive_index / scatterer.wavenumber)
    , m_refractive_index(scatterer.refractive_index)
    , m_medium_index(scatterer.medium_refractive_index)
    , m_wavenumber(scatterer.wavenumber)
{
    validate_scatterer_compatibility();
}

std::map<std::string, std::vector<complex128>> NearFieldEngine::compute_fields(
    const std::vector<std::array<double, 3>>& points,
    const std::vector<std::string>& field_components,
    bool normalize_fields
) const {
    std::map<std::string, std::vector<complex128>> results;

    // Initialize result vectors
    for (const auto& component : field_components) {
        results[component].reserve(points.size());
    }

    // Compute fields at each point
    for (const auto& point : points) {
        double x = point[0], y = point[1], z = point[2];

        // Compute electric and magnetic fields
        auto [Ex, Ey, Ez] = compute_electric_field(x, y, z);
        auto [Hx, Hy, Hz] = compute_magnetic_field(x, y, z);

        // Apply normalization if requested
        if (normalize_fields) {
            auto [E0x, E0y, E0z] = compute_incident_field(x, y, z, "electric");
            auto [H0x, H0y, H0z] = compute_incident_field(x, y, z, "magnetic");

            double E0_mag = std::sqrt(std::norm(E0x) + std::norm(E0y) + std::norm(E0z));
            double H0_mag = std::sqrt(std::norm(H0x) + std::norm(H0y) + std::norm(H0z));

            if (E0_mag > 0) {
                Ex /= E0_mag; Ey /= E0_mag; Ez /= E0_mag;
            }
            if (H0_mag > 0) {
                Hx /= H0_mag; Hy /= H0_mag; Hz /= H0_mag;
            }
        }

        // Store requested field components
        for (const auto& component : field_components) {
            if (component == "Ex") results[component].push_back(Ex);
            else if (component == "Ey") results[component].push_back(Ey);
            else if (component == "Ez") results[component].push_back(Ez);
            else if (component == "Hx") results[component].push_back(Hx);
            else if (component == "Hy") results[component].push_back(Hy);
            else if (component == "Hz") results[component].push_back(Hz);
            else if (component == "|E|") {
                double E_mag = std::sqrt(std::norm(Ex) + std::norm(Ey) + std::norm(Ez));
                results[component].push_back(complex128(E_mag, 0.0));
            }
            else if (component == "|H|") {
                double H_mag = std::sqrt(std::norm(Hx) + std::norm(Hy) + std::norm(Hz));
                results[component].push_back(complex128(H_mag, 0.0));
            }
            else if (component == "|E|²") {
                double E_mag_sq = std::norm(Ex) + std::norm(Ey) + std::norm(Ez);
                results[component].push_back(complex128(E_mag_sq, 0.0));
            }
            else if (component == "|H|²") {
                double H_mag_sq = std::norm(Hx) + std::norm(Hy) + std::norm(Hz);
                results[component].push_back(complex128(H_mag_sq, 0.0));
            }
            else {
                throw std::invalid_argument("Unknown field component: " + component);
            }
        }
    }

    return results;
}

std::tuple<complex128, complex128, complex128> NearFieldEngine::compute_electric_field(
    double x, double y, double z
) const {
    auto [r, theta, phi] = cartesian_to_spherical(x, y, z);

    complex128 Er(0, 0), Etheta(0, 0), Ephi(0, 0);

    // Determine if point is inside or outside scatterer
    bool inside_scatterer = is_point_inside_scatterer(x, y, z);

    // Sum over multipole orders
    for (size_t n = 1; n <= m_max_order; ++n) {
        for (int m = -static_cast<int>(n); m <= static_cast<int>(n); ++m) {
            // Choose appropriate coefficients and Bessel functions
            complex128 coeff_a, coeff_b;
            int bessel_type;

            if (inside_scatterer) {
                // Use internal coefficients cn, dn with regular Bessel functions
                coeff_a = m_scatterer.cn[n-1];  // c_n coefficient
                coeff_b = m_scatterer.dn[n-1];  // d_n coefficient
                bessel_type = 1;  // Regular spherical Bessel functions j_n
            } else {
                // Use scattering coefficients an, bn with outgoing Hankel functions
                coeff_a = m_scatterer.an[n-1];  // a_n coefficient
                coeff_b = m_scatterer.bn[n-1];  // b_n coefficient
                bessel_type = 3;  // Spherical Hankel functions h_n^(1)
            }

            // Compute vector spherical harmonics
            auto [M_vec, N_vec] = compute_vector_spherical_harmonics(n, m, r, theta, phi, bessel_type);
            auto [Mr, Mtheta, Mphi] = M_vec;
            auto [Nr, Ntheta, Nphi] = N_vec;

            // Add multipole contributions to total field
            Er += coeff_a * Mr + coeff_b * Nr;
            Etheta += coeff_a * Mtheta + coeff_b * Ntheta;
            Ephi += coeff_a * Mphi + coeff_b * Nphi;
        }
    }

    // Transform from spherical to Cartesian coordinates
    return spherical_to_cartesian_vector(Er, Etheta, Ephi, theta, phi);
}

std::tuple<complex128, complex128, complex128> NearFieldEngine::compute_magnetic_field(
    double x, double y, double z
) const {
    // Magnetic field computation follows similar pattern to electric field
    // but with different coefficient combinations and normalization
    auto [r, theta, phi] = cartesian_to_spherical(x, y, z);

    complex128 Hr(0, 0), Htheta(0, 0), Hphi(0, 0);

    bool inside_scatterer = is_point_inside_scatterer(x, y, z);
    complex128 impedance = inside_scatterer ?
        complex128(377.0) / m_refractive_index :  // Inside: Z0/n
        complex128(377.0);                        // Outside: Z0

    // Sum over multipole orders (magnetic field has swapped M/N contributions)
    for (size_t n = 1; n <= m_max_order; ++n) {
        for (int m = -static_cast<int>(n); m <= static_cast<int>(n); ++m) {
            complex128 coeff_a, coeff_b;
            int bessel_type;

            if (inside_scatterer) {
                coeff_a = m_scatterer.dn[n-1];  // d_n for magnetic field
                coeff_b = m_scatterer.cn[n-1];  // c_n for magnetic field
                bessel_type = 1;
            } else {
                coeff_a = m_scatterer.bn[n-1];  // b_n for magnetic field
                coeff_b = m_scatterer.an[n-1];  // a_n for magnetic field
                bessel_type = 3;
            }

            auto [M_vec, N_vec] = compute_vector_spherical_harmonics(n, m, r, theta, phi, bessel_type);
            auto [Mr, Mtheta, Mphi] = M_vec;
            auto [Nr, Ntheta, Nphi] = N_vec;

            // Magnetic field has M and N swapped compared to electric field
            Hr += coeff_a * Mr + coeff_b * Nr;
            Htheta += coeff_a * Mtheta + coeff_b * Ntheta;
            Hphi += coeff_a * Mphi + coeff_b * Nphi;
        }
    }

    // Apply impedance normalization and transform coordinates
    Hr /= impedance;
    Htheta /= impedance;
    Hphi /= impedance;

    return spherical_to_cartesian_vector(Hr, Htheta, Hphi, theta, phi);
}

bool NearFieldEngine::is_point_inside_scatterer(double x, double y, double z) const {
    double r = std::sqrt(x*x + y*y + z*z);
    return r < m_radius;
}

double NearFieldEngine::get_field_enhancement(double x, double y, double z, const std::string& field_type) const {
    if (field_type == "electric") {
        auto [Ex, Ey, Ez] = compute_electric_field(x, y, z);
        auto [E0x, E0y, E0z] = compute_incident_field(x, y, z, "electric");

        double E_total = std::sqrt(std::norm(Ex) + std::norm(Ey) + std::norm(Ez));
        double E0_total = std::sqrt(std::norm(E0x) + std::norm(E0y) + std::norm(E0z));

        return E0_total > 0 ? E_total / E0_total : 0.0;
    }
    else if (field_type == "magnetic") {
        auto [Hx, Hy, Hz] = compute_magnetic_field(x, y, z);
        auto [H0x, H0y, H0z] = compute_incident_field(x, y, z, "magnetic");

        double H_total = std::sqrt(std::norm(Hx) + std::norm(Hy) + std::norm(Hz));
        double H0_total = std::sqrt(std::norm(H0x) + std::norm(H0y) + std::norm(H0z));

        return H0_total > 0 ? H_total / H0_total : 0.0;
    }
    else {
        throw std::invalid_argument("field_type must be 'electric' or 'magnetic'");
    }
}

std::tuple<double, double, double> NearFieldEngine::cartesian_to_spherical(double x, double y, double z) const {
    double r = std::sqrt(x*x + y*y + z*z);
    double theta = (r > 0) ? std::acos(z / r) : 0.0;
    double phi = std::atan2(y, x);
    return std::make_tuple(r, theta, phi);
}

std::tuple<complex128, complex128, complex128> NearFieldEngine::spherical_to_cartesian_vector(
    complex128 v_r, complex128 v_theta, complex128 v_phi,
    double theta, double phi
) const {
    double sin_theta = std::sin(theta);
    double cos_theta = std::cos(theta);
    double sin_phi = std::sin(phi);
    double cos_phi = std::cos(phi);

    complex128 v_x = v_r * sin_theta * cos_phi + v_theta * cos_theta * cos_phi - v_phi * sin_phi;
    complex128 v_y = v_r * sin_theta * sin_phi + v_theta * cos_theta * sin_phi + v_phi * cos_phi;
    complex128 v_z = v_r * cos_theta - v_theta * sin_theta;

    return std::make_tuple(v_x, v_y, v_z);
}

std::tuple<complex128, complex128, complex128> NearFieldEngine::compute_incident_field(
    double x, double y, double z, const std::string& field_type
) const {
    // For plane wave incident field (simple case)
    // This should be extended based on the actual source type
    if (field_type == "electric") {
        // Assume plane wave propagating in +z direction, polarized in x
        return std::make_tuple(complex128(1.0, 0.0), complex128(0.0, 0.0), complex128(0.0, 0.0));
    }
    else if (field_type == "magnetic") {
        // Corresponding magnetic field for plane wave
        return std::make_tuple(complex128(0.0, 0.0), complex128(1.0/377.0, 0.0), complex128(0.0, 0.0));
    }
    else {
        throw std::invalid_argument("field_type must be 'electric' or 'magnetic'");
    }
}

void NearFieldEngine::validate_scatterer_compatibility() const {
    // Check if scatterer has cn/dn coefficients computed
    if (m_scatterer.cn.empty() || m_scatterer.dn.empty()) {
        throw std::runtime_error(
            "Near-field computation requires internal field coefficients (cn, dn). "
            "Ensure scatterer is constructed with compute_cn_dn=True."
        );
    }

    // Additional validation can be added here for specific scatterer types
}

// Placeholder for vector spherical harmonics computation
// This is a complex mathematical function that needs careful implementation
std::tuple<
    std::tuple<complex128, complex128, complex128>,
    std::tuple<complex128, complex128, complex128>
> NearFieldEngine::compute_vector_spherical_harmonics(
    int n, int m, double r, double theta, double phi, int bessel_type
) const {
    // This is a placeholder - actual implementation requires:
    // 1. Spherical Bessel functions j_n(kr) or h_n^(1)(kr)
    // 2. Associated Legendre polynomials P_n^m(cos θ)
    // 3. Proper normalization and phase conventions
    // 4. Vector field computation in spherical coordinates

    // TODO: Implement full vector spherical harmonics computation
    // For now, return zero fields to allow compilation
    complex128 zero(0.0, 0.0);
    auto M_vec = std::make_tuple(zero, zero, zero);
    auto N_vec = std::make_tuple(zero, zero, zero);

    return std::make_tuple(M_vec, N_vec);
}
