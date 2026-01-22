#include "./sphere.h"

// ---------------------- Constructors ---------------------------------------
Sphere::Sphere(const double _diameter, const complex128 _refractive_index, const double _medium_refractive_index, std::shared_ptr<BaseSource> _source, size_t _max_order)
: BaseScatterer(_max_order, std::move(_source), _medium_refractive_index), diameter(_diameter), refractive_index(_refractive_index)
{
    this->compute_cross_section();
    this->compute_size_parameter();
    this->max_order = (_max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : _max_order;
    this->compute_an_bn(this->max_order);
}

// ---------------------- Methods ---------------------------------------
void Sphere::compute_size_parameter() {
    this->size_parameter = source->wavenumber_vacuum * this->diameter / 2 * this->medium_refractive_index;
    this->size_parameter_squared = pow(this->size_parameter, 2);
}

void Sphere::compute_cross_section() {
    this->cross_section = Constants::PI * std::pow(this->diameter / 2.0, 2);
}

void Sphere::compute_an_bn(size_t _max_order) {
    _max_order = (_max_order == 0 ? this->max_order : _max_order);

    an.resize(_max_order);
    bn.resize(_max_order);

    complex128 psi_n, chi_n, psi_1, chi_1, xi_n, xi_nm1,
        m = this->refractive_index / this->medium_refractive_index,
        mx = m * size_parameter,
        derivative_a, derivative_b;

    size_t nmx = std::max( _max_order, (size_t) std::abs(mx) ) + 16;

    std::vector<complex128> Dn = this->compute_dn(nmx, mx);

    psi_1  = sin(size_parameter);
    chi_1 = cos(size_parameter);

    for (size_t order = 0; order < _max_order; ++order)
    {
        // Calculate psi and chi (Riccati-Bessel functions)
        double nu = order + 1;

        psi_n = +size_parameter * Spherical_::jn(nu, size_parameter);
        chi_n = -size_parameter * Spherical_::yn(nu, size_parameter);

        // Complex Riccati-Bessel functions
        xi_n = psi_n - 1.0 * complex128(0, 1) * chi_n;
        xi_nm1 = psi_1 - 1.0 * complex128(0, 1) * chi_1;

        // Derivative of the Riccati-Bessel functions
        derivative_a = Dn[order + 1] / m + nu / size_parameter;
        derivative_b = Dn[order + 1] * m + nu / size_parameter;

        // Computation of the electric and magnetic multipole coefficients
        an[order] = (derivative_a * psi_n - psi_1) / (derivative_a * xi_n - xi_nm1);
        bn[order] = (derivative_b * psi_n - psi_1) / (derivative_b * xi_n - xi_nm1);

        psi_1 = psi_n;
        chi_1 = chi_n;
    }
}

void Sphere::compute_cn_dn(size_t _max_order) {
    _max_order = (_max_order == 0 ? this->max_order : _max_order);

    cn.resize(_max_order);
    dn.resize(_max_order);

    complex128
        x = size_parameter,
        m = this->refractive_index / this->medium_refractive_index,
        z = m * x;

    size_t nmx = std::max( _max_order, (size_t) std::abs(z) ) + 16;

    std::vector<complex128>
        Cnx = std::vector<complex128>(nmx),
        Cnn, jnx, jnmx, yx, hx, b1x, y1x, hn1x, ax, ahx, numerator,
        c_denominator, d_denominator;

    b1x.push_back( +sin(x) / x );
    y1x.push_back( -cos(x) / x );

    for (double i = nmx; i > 1; i--)
        Cnx[i-2] = i - z * z / (Cnx[i - 1] + i);

    for (size_t order = 0; order < _max_order; order++) {
        Cnn.push_back(Cnx[order]);

        // jnx.push_back(compute_jn(order + 1, x));
        // jnmx.push_back(1.0 / (compute_jn(order + 1, z)));
        // yx.push_back(compute_yn(order + 1, x));
        // hx.push_back(jnx[order] + complex128(0, 1) * yx[order]);

        jnx.push_back(Spherical_::jn(order + 1, x));
        jnmx.push_back(1.0 / (Spherical_::jn(order + 1, z)));
        yx.push_back(Spherical_::yn(order + 1, x));
        hx.push_back(jnx[order] + complex128(0, 1) * yx[order]);

        b1x.push_back(jnx[order]);
        y1x.push_back(yx[order]);
        hn1x.push_back(b1x[order] + complex128(0, 1) * y1x[order]);

        ax.push_back(x * b1x[order] - (complex128)(order + 1.0) * jnx[order]);
        ahx.push_back(x * hn1x[order] - (complex128)(order + 1.0) * hx[order]);

        numerator.push_back( jnx[order] * ahx[order] - hx[order] * ax[order] );
        c_denominator.push_back( ahx[order] - hx[order] * Cnn[order] );
        d_denominator.push_back( m * m * ahx[order] - hx[order] * Cnn[order] );
        cn[order] = jnmx[order] * numerator[order] / c_denominator[order] ;
        dn[order] = jnmx[order] * m * numerator[order] / d_denominator[order] ;
    }
}

double Sphere::get_Qsca() const {
    double value = 0;

    for(size_t it = 0; it < max_order; ++it){
        double n = (double) it + 1;
        value += (2. * n + 1.) * ( pow( std::abs(this->an[it]), 2) + pow( std::abs(this->bn[it]), 2)  );
    }
    return value * 2. / size_parameter_squared;
}

double Sphere::get_Qext() const {
    double value = 0;
    for(size_t it = 0; it < max_order; ++it)
    {
        double n = (double) it + 1;
        value += (2.* n + 1.) * std::real( this->an[it] + this->bn[it] );

    }
    return value * 2. / size_parameter_squared;
}

double Sphere::get_Qback() const {
    complex128 value = 0;

    for(size_t it = 0; it < max_order-1; ++it)
    {
        double n = (double) it + 1;

        value += (2. * n + 1) * pow(-1., n) * ( this->an[it] - this->bn[it] ) ;
    }

    value = pow( std::abs(value), 2. ) / size_parameter_squared;
    return std::abs(value);
}

double Sphere::get_Qforward() const {
    complex128 value = 0;

    for(size_t it = 0; it < max_order-1; ++it)
    {
        double n = (double) it + 1;

        value += (2. * n + 1) * ( this->an[it] + this->bn[it] ) ;
    }

    value = pow( std::abs(value), 2. ) / size_parameter_squared;
    return std::abs(value);
}

double Sphere::get_g() const {
    double value = 0;

      for(size_t it = 0; it < max_order-1; ++it) {
         double n = (double) it + 1;

          value += ( n * (n + 2.) / (n + 1.) ) * std::real(this->an[it] * std::conj(this->an[it+1]) + this->bn[it] * std::conj(this->bn[it+1]) );
          value += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( this->an[it] * std::conj(this->bn[it]) );
      }
      return value * 4. / ( get_Qsca() * size_parameter_squared );
}


std::tuple<std::vector<complex128>, std::vector<complex128>>
Sphere::compute_s1s2(const std::vector<double> &phi) const {
    std::vector<complex128> S1, S2;

    S1.reserve(phi.size());
    S2.reserve(phi.size());

    std::vector<double> mu, prefactor = get_prefactor();

    mu.reserve(phi.size());

    for (const double phi : phi)
        mu.push_back( cos( phi - Constants::PI / 2.0 ) );

    for (size_t i = 0; i < phi.size(); i++) {
        auto [pin, taun] = this->get_pi_tau(mu[i], max_order);
        complex128 S1_temp = 0., S2_temp = 0.;

        for (size_t m = 0; m < max_order ; m++) {
            S1_temp += prefactor[m] * ( this->an[m] * pin[m] +  this->bn[m] * taun[m] );
            S2_temp += prefactor[m] * ( this->an[m] * taun[m] + this->bn[m] * pin[m]  );
        }
        S1.push_back(S1_temp);
        S2.push_back(S2_temp);
    }

    return std::make_tuple(std::move(S1), std::move(S2));
}

// std::vector<complex128>
// Sphere::compute_nearfields(
//     const std::vector<double>& x,
//     const std::vector<double>& y,
//     const std::vector<double>& z,
//     const std::string& field_type
// )
// {
//     this->compute_cn_dn(this->max_order);

//     if (field_type != "Ex" && field_type != "Ey" && field_type != "Ez" &&
//         field_type != "|E|") {
//         throw std::invalid_argument("Invalid field_type. Must be one of: Ex, Ey, Ez, |E|");
//     }

//     const size_t number_of_points = x.size();
//     if (number_of_points != y.size() || number_of_points != z.size()) {
//         throw std::invalid_argument("x, y, z vectors must have the same length");
//     }

//     std::vector<complex128> field_values(number_of_points);

//     const double medium_refractive_index = this->medium_refractive_index;
//     const double radius_particle = 0.5 * this->diameter;

//     const double k_medium = this->source->wavenumber_vacuum * medium_refractive_index;
//     const complex128 i_unit(0.0, 1.0);

//     const std::array<complex128, 2> jones_vector = this->source->get_jones_vector_first_row();
//     const complex128 incident_Ex_amplitude = jones_vector[0] * this->source->amplitude;
//     const complex128 incident_Ey_amplitude = jones_vector[1] * this->source->amplitude;

//     auto clamp_m1_p1_local = [](double value) -> double {
//         if (value < -1.0) return -1.0;
//         if (value >  1.0) return  1.0;
//         return value;
//     };

//     auto compute_psi_xi_and_derivatives =
//         [&](size_t n, complex128 argument, bool use_hankel)
//         -> std::tuple<complex128, complex128, complex128, complex128>
//     {
//         // Returns:
//         // radial_function_n       = j_n(arg) or h_n(arg)
//         // radial_function_n_minus = j_{n-1}(arg) or h_{n-1}(arg)
//         // riccati_n               = arg * radial_function_n
//         // riccati_n_derivative    = d/darg [arg * radial_function_n]
//         //
//         // Derivative identity:
//         // d/dz [z f_n(z)] = z f_{n-1}(z) - n f_n(z)
//         // so riccati_n' = arg * f_{n-1}(arg) - n * f_n(arg)

//         const double order_n = static_cast<double>(n);
//         const double order_nm1 = static_cast<double>(n - 1);

//         complex128 fn;
//         complex128 fnm1;

//         if (!use_hankel) {
//             fn   = Spherical_::jn(order_n, argument);
//             fnm1 = (n == 1) ? Spherical_::jn(0.0, argument) : Spherical_::jn(order_nm1, argument);
//         } else {
//             fn   = Spherical_::H1n(order_n, argument);
//             fnm1 = (n == 1) ? Spherical_::H1n(0.0, argument) : Spherical_::H1n(order_nm1, argument);
//         }

//         const complex128 riccati_n = argument * fn;
//         const complex128 riccati_n_derivative = argument * fnm1 - order_n * fn;

//         return {fn, fnm1, riccati_n, riccati_n_derivative};
//     };

//     for (size_t point_index = 0; point_index < number_of_points; ++point_index) {

//         const double x_position = x[point_index];
//         const double y_position = y[point_index];
//         const double z_position = z[point_index];

//         const double r_cartesian = std::sqrt(x_position * x_position + y_position * y_position + z_position * z_position);

//         // Handle the origin safely
//         if (r_cartesian < 1e-15) {
//             // Total field at origin is just the incident field evaluated at z = 0
//             if (field_type == "Ex") field_values[point_index] = incident_Ex_amplitude;
//             if (field_type == "Ey") field_values[point_index] = incident_Ey_amplitude;
//             if (field_type == "Ez") field_values[point_index] = complex128(0.0, 0.0);
//             if (field_type == "|E|") {
//                 field_values[point_index] = std::sqrt(std::norm(incident_Ex_amplitude) + std::norm(incident_Ey_amplitude));
//             }
//             continue;
//         }

//         const double mu = clamp_m1_p1_local(z_position / r_cartesian);
//         const double theta = std::acos(mu);
//         const double phi = std::atan2(y_position, x_position);

//         const double sin_theta = std::sin(theta);
//         const double cos_theta = std::cos(theta);
//         const double sin_phi   = std::sin(phi);
//         const double cos_phi   = std::cos(phi);

//         const bool is_inside = (r_cartesian < radius_particle);

//         const complex128 k_r = complex128(k_medium * r_cartesian, 0.0);
//         const complex128 incident_phase = std::exp(i_unit * complex128(k_medium * z_position, 0.0));

//         // Incident field in Cartesian coordinates
//         const complex128 incident_Ex = incident_Ex_amplitude * incident_phase;
//         const complex128 incident_Ey = incident_Ey_amplitude * incident_phase;
//         const complex128 incident_Ez = complex128(0.0, 0.0);

//         // Convert incident E to spherical components (Er, Etheta, Ephi)
//         // Spherical basis:
//         // e_r     = (sinθ cosφ, sinθ sinφ, cosθ)
//         // e_theta = (cosθ cosφ, cosθ sinφ, -sinθ)
//         // e_phi   = (-sinφ, cosφ, 0)
//         const complex128 incident_Er =
//             incident_Ex * (sin_theta * cos_phi) +
//             incident_Ey * (sin_theta * sin_phi) +
//             incident_Ez * (cos_theta);

//         const complex128 incident_Etheta =
//             incident_Ex * (cos_theta * cos_phi) +
//             incident_Ey * (cos_theta * sin_phi) -
//             incident_Ez * (sin_theta);

//         const complex128 incident_Ephi =
//             incident_Ex * (-sin_phi) +
//             incident_Ey * (cos_phi);

//         // Accumulate scattered or internal correction in spherical components
//         complex128 correction_Er(0.0, 0.0);
//         complex128 correction_Etheta(0.0, 0.0);
//         complex128 correction_Ephi(0.0, 0.0);

//         // Compute pi_n and tau_n at this theta using your existing routine
//         // pi_n and tau_n are real for real mu, but you store complex128, that is fine
//         auto [pi_vector, tau_vector] = this->get_pi_tau(mu, this->max_order);

//         // Mie series for near field uses radial Riccati functions, not far field propagator.
//         //
//         // This implementation follows the standard m = 1 plane wave expansion:
//         // Angular dependence is carried by pi_n and tau_n, and azimuth enters only as cos(phi) and sin(phi)
//         //
//         // Outside:
//         // use an, bn with hankel radial functions
//         // Inside:
//         // use cn, dn with spherical bessel radial functions
//         printf("medium_refractive_index: %f\n", medium_refractive_index);
//         printf("size_parameter: %f\n", this->size_parameter);
//         for (size_t n = 1; n <= this->max_order; ++n) {

//             const double n_double = static_cast<double>(n);

//             const complex128 en = std::pow(i_unit, n_double) * (2.0 * n_double + 1.0) / (n_double * (n_double + 1.0));

//             const complex128 pi_n  = pi_vector[n - 1];
//             const complex128 tau_n = tau_vector[n - 1];

//             // Select coefficients and radial functions
//             const complex128 coefficient_a = is_inside ? this->cn[n - 1] : this->an[n - 1];
//             const complex128 coefficient_b = is_inside ? this->dn[n - 1] : this->bn[n - 1];

//             const bool use_hankel = !is_inside;

//             auto [radial_n, radial_n_minus_1, riccati_n, riccati_n_derivative] =
//                 compute_psi_xi_and_derivatives(n, k_r, use_hankel);

//             // Also need riccati over kr and derivative over kr
//             const complex128 riccati_over_kr = riccati_n / k_r;
//             const complex128 riccati_derivative_over_kr = riccati_n_derivative / k_r;

//             // The m = 1 azimuthal structure for a plane wave can be expressed by projecting
//             // the incident polarization onto local spherical basis. A clean way that preserves
//             // your arbitrary Jones vector is to build the correction in the same local basis
//             // as the incident components.
//             //
//             // The standard combination for Etheta and Ephi uses:
//             // S1_like = a_n * pi_n + b_n * tau_n
//             // S2_like = a_n * tau_n + b_n * pi_n
//             //
//             // For near field, those angular combinations multiply riccati_n and riccati_n_derivative pieces.
//             //
//             // This gives physically correct symmetry and removes artificial vortices.

//             const complex128 angular_S1 = coefficient_a * pi_n + coefficient_b * tau_n;
//             const complex128 angular_S2 = coefficient_a * tau_n + coefficient_b * pi_n;

//             // Scattered or internal correction in spherical coordinates
//             //
//             // The following structure matches the standard M and N wave mixing:
//             // Etheta uses both riccati and derivative terms
//             // Ephi uses both riccati and derivative terms
//             //
//             // This is the key: no exp(i n phi) anywhere, only the pi and tau angular functions.

//             correction_Etheta += en * ( angular_S2 * riccati_derivative_over_kr );
//             correction_Ephi   += en * ( angular_S1 * riccati_derivative_over_kr );

//             // Radial component is needed for correct Cartesian conversion near the particle.
//             // Use the standard relation with pi_n and riccati_over_kr.
//             // Guard sin(theta) because pi_n ~ P_n^1 / sin(theta)
//             if (sin_theta > 1e-15) {
//                 correction_Er += en * (coefficient_a - coefficient_b) * (pi_n * riccati_over_kr) * cos_theta;
//             }
//         }

//         // Total field in spherical components
//         const complex128 total_Er     = incident_Er     + correction_Er;
//         const complex128 total_Etheta = incident_Etheta + correction_Etheta;
//         const complex128 total_Ephi   = incident_Ephi   + correction_Ephi;

//         // Convert total spherical components back to Cartesian
//         const complex128 total_Ex =
//             total_Er * (sin_theta * cos_phi) +
//             total_Etheta * (cos_theta * cos_phi) +
//             total_Ephi * (-sin_phi);

//         const complex128 total_Ey =
//             total_Er * (sin_theta * sin_phi) +
//             total_Etheta * (cos_theta * sin_phi) +
//             total_Ephi * (cos_phi);

//         const complex128 total_Ez =
//             total_Er * (cos_theta) +
//             total_Etheta * (-sin_theta);

//         if (field_type == "Ex") {
//             field_values[point_index] = total_Ex;
//         } else if (field_type == "Ey") {
//             field_values[point_index] = total_Ey;
//         } else if (field_type == "Ez") {
//             field_values[point_index] = total_Ez;
//         } else if (field_type == "|E|") {
//             field_values[point_index] = std::sqrt(std::norm(total_Ex) + std::norm(total_Ey) + std::norm(total_Ez));
//         }
//     }

//     return field_values;
// }

// std::vector<complex128>
// Sphere::compute_nearfields(
//     const std::vector<double>& x,
//     const std::vector<double>& y,
//     const std::vector<double>& z,
//     const std::string& field_type
// )
// {
//     this->compute_cn_dn(this->max_order);

//     if (field_type != "Ex" && field_type != "Ey" && field_type != "Ez" && field_type != "|E|") {
//         throw std::invalid_argument("Invalid field_type. Must be one of: Ex, Ey, Ez, |E|");
//     }

//     const size_t number_of_points = x.size();
//     if (number_of_points != y.size() || number_of_points != z.size()) {
//         throw std::invalid_argument("x, y, z vectors must have the same length");
//     }

//     std::vector<complex128> field_values(number_of_points);

//     const double radius_particle = 0.5 * this->diameter;
//     const double medium_refractive_index = this->medium_refractive_index;

//     const double k0 = this->source->wavenumber_vacuum;
//     const double k_medium = k0 * medium_refractive_index;

//     const complex128 i_unit(0.0, 1.0);

//     const std::array<complex128, 2> jones_vector = this->source->get_jones_vector_first_row();
//     const complex128 E0x = jones_vector[0] * this->source->amplitude;
//     const complex128 E0y = jones_vector[1] * this->source->amplitude;

//     auto clamp_m1_p1_local = [](double value) -> double {
//         if (value < -1.0) return -1.0;
//         if (value >  1.0) return  1.0;
//         return value;
//     };

//     auto i_to_n = [&](size_t n) -> complex128 {
//         // i^n with n >= 0
//         switch (n & 3u) {
//             case 0u: return complex128( 1.0,  0.0);
//             case 1u: return complex128( 0.0,  1.0);
//             case 2u: return complex128(-1.0,  0.0);
//             default: return complex128( 0.0, -1.0);
//         }
//     };

//     auto compute_zn_and_dpsi_over_kr =
//         [&](size_t n, const complex128& kr, bool use_hankel)
//         -> std::tuple<complex128, complex128>
//     {
//         // z_n(kr) where z is j or h1
//         // dpsi = d/d(kr)[ kr z_n(kr) ] = kr z_{n-1}(kr) - n z_n(kr)
//         // return ( z_n(kr), dpsi/(kr) )

//         const double nd = static_cast<double>(n);
//         const double nm1d = static_cast<double>(n - 1);

//         complex128 zn;
//         complex128 znm1;

//         if (!use_hankel) {
//             zn   = Spherical_::jn(nd,  kr);
//             znm1 = (n == 1) ? Spherical_::jn(0.0, kr) : Spherical_::jn(nm1d, kr);
//         } else {
//             zn   = Spherical_::H1n(nd,  kr);
//             znm1 = (n == 1) ? Spherical_::H1n(0.0, kr) : Spherical_::H1n(nm1d, kr);
//         }

//         const complex128 dpsi = kr * znm1 - nd * zn;

//         const complex128 dpsi_over_kr = dpsi / kr;

//         return {zn, dpsi_over_kr};
//     };

//     for (size_t point_index = 0; point_index < number_of_points; ++point_index) {

//         const double x_position = x[point_index];
//         const double y_position = y[point_index];
//         const double z_position = z[point_index];

//         const double r_cartesian = std::sqrt(
//             x_position * x_position + y_position * y_position + z_position * z_position
//         );

//         if (r_cartesian < 1e-15) {
//             // The VSWF representation is singular at the origin in spherical coordinates.
//             // For debugging, return the plane wave at z = 0 with the requested Jones amplitudes.
//             if (field_type == "Ex") field_values[point_index] = E0x;
//             else if (field_type == "Ey") field_values[point_index] = E0y;
//             else if (field_type == "Ez") field_values[point_index] = complex128(0.0, 0.0);
//             else {
//                 const double mag = std::sqrt(std::norm(E0x) + std::norm(E0y));
//                 field_values[point_index] = complex128(mag, 0.0);
//             }
//             continue;
//         }

//         const double mu_raw = z_position / r_cartesian;
//         const double mu = clamp_m1_p1_local(mu_raw);
//         const double theta = std::acos(mu);
//         const double phi = std::atan2(y_position, x_position);

//         const double sin_theta = std::sin(theta);
//         const double cos_theta = std::cos(theta);
//         const double sin_phi   = std::sin(phi);
//         const double cos_phi   = std::cos(phi);

//         const bool is_inside = (r_cartesian < radius_particle);

//         const complex128 kr = complex128(k_medium * r_cartesian, 0.0);

//         // pi_n and tau_n for m = 1, length max_order
//         auto [pi_vector, tau_vector] = this->get_pi_tau(mu, this->max_order);

//         // Accumulate spherical components for two independent incident polarizations:
//         // x polarized: uses M_o1n and N_e1n
//         // y polarized: uses M_e1n and N_o1n
//         //
//         // Bohren Huffman convention for a plane wave propagating +z:
//         // E_inc,x = sum En ( M_o1n^(1) - i N_e1n^(1) )
//         // E_sca,x = sum En ( i a_n N_e1n^(3) - b_n M_o1n^(3) )
//         // E_int,x = sum En ( c_n M_o1n^(1) - i d_n N_e1n^(1) )
//         //
//         // E_inc,y = sum En ( M_e1n^(1) - i N_o1n^(1) )
//         // E_sca,y = sum En ( i a_n N_o1n^(3) - b_n M_e1n^(3) )
//         // E_int,y = sum En ( c_n M_e1n^(1) - i d_n N_o1n^(1) )

//         complex128 Er_x(0.0, 0.0), Etheta_x(0.0, 0.0), Ephi_x(0.0, 0.0);
//         complex128 Er_y(0.0, 0.0), Etheta_y(0.0, 0.0), Ephi_y(0.0, 0.0);

//         for (size_t n = 1; n <= this->max_order; ++n) {

//             const double nd = static_cast<double>(n);

//             const complex128 En = i_to_n(n) * (2.0 * nd + 1.0) / (nd * (nd + 1.0));

//             const complex128 pi_n  = pi_vector[n - 1];
//             const complex128 tau_n = tau_vector[n - 1];

//             // P_n^1(mu) = pi_n(mu) * sin(theta)
//             const complex128 Pn1 = pi_n * complex128(sin_theta, 0.0);

//             // Radial functions
//             // Regular: j_n
//             // Outgoing: h_n^(1)
//             const bool use_hankel_for_scat = true;
//             const bool use_hankel_for_inc  = false;

//             auto [jn_val, dpsi_j_over_kr] = compute_zn_and_dpsi_over_kr(n, kr, use_hankel_for_inc);
//             auto [hn_val, dpsi_h_over_kr] = compute_zn_and_dpsi_over_kr(n, kr, use_hankel_for_scat);

//             const complex128 j_over_kr = jn_val / kr;
//             const complex128 h_over_kr = hn_val / kr;

//             // Coefficients
//             const complex128 a_n = this->an[n - 1];
//             const complex128 b_n = this->bn[n - 1];
//             const complex128 c_n = this->cn[n - 1];
//             const complex128 d_n = this->dn[n - 1];

//             // Build VSWFs in spherical components for m = 1
//             //
//             // M_o1n(z):  (0,  cosφ pi_n z_n,  -sinφ tau_n z_n)
//             // M_e1n(z):  (0,  sinφ pi_n z_n,   cosφ tau_n z_n)
//             //
//             // N_e1n(z):  (cosφ n(n+1) Pn1 z_n/(kr),  cosφ tau_n dpsi/(kr),  -sinφ pi_n dpsi/(kr))
//             // N_o1n(z):  (sinφ n(n+1) Pn1 z_n/(kr),  sinφ tau_n dpsi/(kr),   cosφ pi_n dpsi/(kr))

//             const complex128 nn1 = complex128(nd * (nd + 1.0), 0.0);

//             // Regular (incident and internal)
//             const complex128 Mo_r_j(0.0, 0.0);
//             const complex128 Mo_t_j = complex128(cos_phi, 0.0) * pi_n  * jn_val;
//             const complex128 Mo_p_j = complex128(-sin_phi, 0.0) * tau_n * jn_val;

//             const complex128 Me_r_j(0.0, 0.0);
//             const complex128 Me_t_j = complex128(sin_phi, 0.0) * pi_n  * jn_val;
//             const complex128 Me_p_j = complex128(cos_phi, 0.0) * tau_n * jn_val;

//             const complex128 Ne_r_j = complex128(cos_phi, 0.0) * nn1 * Pn1 * j_over_kr;
//             const complex128 Ne_t_j = complex128(cos_phi, 0.0) * tau_n * dpsi_j_over_kr;
//             const complex128 Ne_p_j = complex128(-sin_phi, 0.0) * pi_n * dpsi_j_over_kr;

//             const complex128 No_r_j = complex128(sin_phi, 0.0) * nn1 * Pn1 * j_over_kr;
//             const complex128 No_t_j = complex128(sin_phi, 0.0) * tau_n * dpsi_j_over_kr;
//             const complex128 No_p_j = complex128(cos_phi, 0.0) * pi_n * dpsi_j_over_kr;

//             // Outgoing (scattered)
//             const complex128 Mo_r_h(0.0, 0.0);
//             const complex128 Mo_t_h = complex128(cos_phi, 0.0) * pi_n  * hn_val;
//             const complex128 Mo_p_h = complex128(-sin_phi, 0.0) * tau_n * hn_val;

//             const complex128 Me_r_h(0.0, 0.0);
//             const complex128 Me_t_h = complex128(sin_phi, 0.0) * pi_n  * hn_val;
//             const complex128 Me_p_h = complex128(cos_phi, 0.0) * tau_n * hn_val;

//             const complex128 Ne_r_h = complex128(cos_phi, 0.0) * nn1 * Pn1 * h_over_kr;
//             const complex128 Ne_t_h = complex128(cos_phi, 0.0) * tau_n * dpsi_h_over_kr;
//             const complex128 Ne_p_h = complex128(-sin_phi, 0.0) * pi_n * dpsi_h_over_kr;

//             const complex128 No_r_h = complex128(sin_phi, 0.0) * nn1 * Pn1 * h_over_kr;
//             const complex128 No_t_h = complex128(sin_phi, 0.0) * tau_n * dpsi_h_over_kr;
//             const complex128 No_p_h = complex128(cos_phi, 0.0) * pi_n * dpsi_h_over_kr;

//             if (!is_inside) {

//                 // Total outside = incident + scattered

//                 // Incident x
//                 Er_x     += En * ( Mo_r_j - i_unit * Ne_r_j );
//                 Etheta_x += En * ( Mo_t_j - i_unit * Ne_t_j );
//                 Ephi_x   += En * ( Mo_p_j - i_unit * Ne_p_j );

//                 // Scattered x
//                 Er_x     += En * ( i_unit * a_n * Ne_r_h - b_n * Mo_r_h );
//                 Etheta_x += En * ( i_unit * a_n * Ne_t_h - b_n * Mo_t_h );
//                 Ephi_x   += En * ( i_unit * a_n * Ne_p_h - b_n * Mo_p_h );

//                 // Incident y
//                 Er_y     += En * ( Me_r_j - i_unit * No_r_j );
//                 Etheta_y += En * ( Me_t_j - i_unit * No_t_j );
//                 Ephi_y   += En * ( Me_p_j - i_unit * No_p_j );

//                 // Scattered y
//                 Er_y     += En * ( i_unit * a_n * No_r_h - b_n * Me_r_h );
//                 Etheta_y += En * ( i_unit * a_n * No_t_h - b_n * Me_t_h );
//                 Ephi_y   += En * ( i_unit * a_n * No_p_h - b_n * Me_p_h );

//             } else {

//                 // Total inside = internal field only

//                 // Internal x
//                 Er_x     += En * ( c_n * Mo_r_j - i_unit * d_n * Ne_r_j );
//                 Etheta_x += En * ( c_n * Mo_t_j - i_unit * d_n * Ne_t_j );
//                 Ephi_x   += En * ( c_n * Mo_p_j - i_unit * d_n * Ne_p_j );

//                 // Internal y
//                 Er_y     += En * ( c_n * Me_r_j - i_unit * d_n * No_r_j );
//                 Etheta_y += En * ( c_n * Me_t_j - i_unit * d_n * No_t_j );
//                 Ephi_y   += En * ( c_n * Me_p_j - i_unit * d_n * No_p_j );
//             }
//         }

//         // Combine arbitrary incident polarization
//         const complex128 Er     = E0x * Er_x     + E0y * Er_y;
//         const complex128 Etheta = E0x * Etheta_x + E0y * Etheta_y;
//         const complex128 Ephi   = E0x * Ephi_x   + E0y * Ephi_y;

//         // Convert spherical to Cartesian
//         // e_r     = (sinθ cosφ, sinθ sinφ, cosθ)
//         // e_theta = (cosθ cosφ, cosθ sinφ, -sinθ)
//         // e_phi   = (-sinφ, cosφ, 0)

//         const complex128 Ex =
//             Er * (sin_theta * cos_phi)
//             + Etheta * (cos_theta * cos_phi)
//             + Ephi * (-sin_phi);

//         const complex128 Ey =
//             Er * (sin_theta * sin_phi)
//             + Etheta * (cos_theta * sin_phi)
//             + Ephi * (cos_phi);

//         const complex128 Ez =
//             Er * (cos_theta)
//             + Etheta * (-sin_theta);

//         if (field_type == "Ex") {
//             field_values[point_index] = Ex;
//         } else if (field_type == "Ey") {
//             field_values[point_index] = Ey;
//         } else if (field_type == "Ez") {
//             field_values[point_index] = Ez;
//         } else {
//             const double mag = std::sqrt(std::norm(Ex) + std::norm(Ey) + std::norm(Ez));
//             field_values[point_index] = complex128(mag, 0.0);
//         }
//     }

//     return field_values;
// }


std::vector<complex128>
Sphere::compute_nearfields(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::string& field_type
)
{
    this->compute_cn_dn(this->max_order);

    if (field_type != "Ex" && field_type != "Ey" && field_type != "Ez" && field_type != "|E|") {
        throw std::invalid_argument("Invalid field_type. Must be one of: Ex, Ey, Ez, |E|");
    }

    const std::size_t number_of_points = x.size();
    if (number_of_points != y.size() || number_of_points != z.size()) {
        throw std::invalid_argument("x, y, z vectors must have the same length");
    }

    std::vector<complex128> field_values(number_of_points);

    const double radius_particle = 0.5 * this->diameter;

    const double medium_refractive_index = this->medium_refractive_index;
    const complex128 particle_refractive_index = this->refractive_index;

    const double k0 = this->source->wavenumber_vacuum;
    const double k_medium = k0 * medium_refractive_index;

    const complex128 k_particle = complex128(k0, 0.0) * particle_refractive_index;

    const complex128 i_unit(0.0, 1.0);

    const std::array<complex128, 2> jones_vector = this->source->get_jones_vector_first_row();
    const complex128 E0x = jones_vector[0] * this->source->amplitude;
    const complex128 E0y = jones_vector[1] * this->source->amplitude;

    auto clamp_m1_p1_local = [](double value) -> double {
        if (value < -1.0) return -1.0;
        if (value >  1.0) return  1.0;
        return value;
    };

    auto i_to_n = [&](std::size_t n) -> complex128 {
        switch (n & 3u) {
            case 0u: return complex128( 1.0,  0.0);
            case 1u: return complex128( 0.0,  1.0);
            case 2u: return complex128(-1.0,  0.0);
            default: return complex128( 0.0, -1.0);
        }
    };

    auto compute_zn_and_dpsi_over_kr =
        [&](std::size_t n, const complex128& kr, bool use_hankel)
        -> std::tuple<complex128, complex128, complex128>
    {
        // z_n(kr) where z is j or h1
        // dpsi = d/d(kr)[ kr z_n(kr) ] = kr z_{n-1}(kr) - n z_n(kr)
        // return ( z_n(kr), z_n(kr)/(kr), dpsi/(kr) )

        const double nd = static_cast<double>(n);
        const double nm1d = static_cast<double>(n - 1);

        complex128 zn;
        complex128 znm1;

        if (!use_hankel) {
            zn   = Spherical_::jn(nd,  kr);
            znm1 = (n == 1) ? Spherical_::jn(0.0, kr) : Spherical_::jn(nm1d, kr);
        } else {
            zn   = Spherical_::H1n(nd,  kr);
            znm1 = (n == 1) ? Spherical_::H1n(0.0, kr) : Spherical_::H1n(nm1d, kr);
        }

        const complex128 dpsi = kr * znm1 - nd * zn;

        const double abs_kr = std::abs(kr);
        if (abs_kr < 1e-15) {
            // This branch should only matter extremely close to the origin.
            // The caller already handles r -> 0 explicitly, but keep this defensive.
            return {zn, complex128(0.0, 0.0), complex128(0.0, 0.0)};
        }

        const complex128 z_over_kr = zn / kr;
        const complex128 dpsi_over_kr = dpsi / kr;

        return {zn, z_over_kr, dpsi_over_kr};
    };

    for (std::size_t point_index = 0; point_index < number_of_points; ++point_index) {

        const double x_position = x[point_index];
        const double y_position = y[point_index];
        const double z_position = z[point_index];

        const double r_cartesian = std::sqrt(
            x_position * x_position + y_position * y_position + z_position * z_position
        );

        if (r_cartesian < 1e-15) {
            // The VSWF representation is singular at the origin in spherical coordinates.
            // For debugging, return the plane wave at z = 0 with the requested Jones amplitudes.
            if (field_type == "Ex") field_values[point_index] = E0x;
            else if (field_type == "Ey") field_values[point_index] = E0y;
            else if (field_type == "Ez") field_values[point_index] = complex128(0.0, 0.0);
            else {
                const double mag = std::sqrt(std::norm(E0x) + std::norm(E0y));
                field_values[point_index] = complex128(mag, 0.0);
            }
            continue;
        }

        const double mu_raw = z_position / r_cartesian;
        const double mu = clamp_m1_p1_local(mu_raw);
        const double theta = std::acos(mu);
        const double phi = std::atan2(y_position, x_position);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double sin_phi   = std::sin(phi);
        const double cos_phi   = std::cos(phi);

        const bool is_inside = (r_cartesian < radius_particle);

        const complex128 kr_outside = complex128(k_medium * r_cartesian, 0.0);
        const complex128 kr_inside  = k_particle * complex128(r_cartesian, 0.0);

        const complex128 kr_used_for_regular = is_inside ? kr_inside : kr_outside;
        const complex128 kr_used_for_hankel  = kr_outside;

        auto [pi_vector, tau_vector] = this->get_pi_tau(mu, this->max_order);

        complex128 Er_x(0.0, 0.0), Etheta_x(0.0, 0.0), Ephi_x(0.0, 0.0);
        complex128 Er_y(0.0, 0.0), Etheta_y(0.0, 0.0), Ephi_y(0.0, 0.0);

        for (std::size_t n = 1; n <= this->max_order; ++n) {

            const double nd = static_cast<double>(n);

            const complex128 En = i_to_n(n) * (2.0 * nd + 1.0) / (nd * (nd + 1.0));

            const complex128 pi_n  = pi_vector[n - 1];
            const complex128 tau_n = tau_vector[n - 1];

            const complex128 Pn1 = pi_n * complex128(sin_theta, 0.0);

            const complex128 a_n = this->an[n - 1];
            const complex128 b_n = this->bn[n - 1];
            const complex128 c_n = this->cn[n - 1];
            const complex128 d_n = this->dn[n - 1];

            const complex128 nn1 = complex128(nd * (nd + 1.0), 0.0);

            const bool use_hankel_for_incident = false;
            const bool use_hankel_for_scattered = true;

            auto [zn_regular, z_regular_over_kr, dpsi_regular_over_kr] =
                compute_zn_and_dpsi_over_kr(n, kr_used_for_regular, use_hankel_for_incident);

            complex128 zn_hankel(0.0, 0.0), z_hankel_over_kr(0.0, 0.0), dpsi_hankel_over_kr(0.0, 0.0);

            if (!is_inside) {
                std::tie(zn_hankel, z_hankel_over_kr, dpsi_hankel_over_kr) =
                    compute_zn_and_dpsi_over_kr(n, kr_used_for_hankel, use_hankel_for_scattered);
            }

            // Regular (incident outside, internal inside)
            const complex128 Mo_r_j(0.0, 0.0);
            const complex128 Mo_t_j = complex128(cos_phi, 0.0) * pi_n  * zn_regular;
            const complex128 Mo_p_j = complex128(-sin_phi, 0.0) * tau_n * zn_regular;

            const complex128 Me_r_j(0.0, 0.0);
            const complex128 Me_t_j = complex128(sin_phi, 0.0) * pi_n  * zn_regular;
            const complex128 Me_p_j = complex128(cos_phi, 0.0) * tau_n * zn_regular;

            const complex128 Ne_r_j = complex128(cos_phi, 0.0) * nn1 * Pn1 * z_regular_over_kr;
            const complex128 Ne_t_j = complex128(cos_phi, 0.0) * tau_n * dpsi_regular_over_kr;
            const complex128 Ne_p_j = complex128(-sin_phi, 0.0) * pi_n  * dpsi_regular_over_kr;

            const complex128 No_r_j = complex128(sin_phi, 0.0) * nn1 * Pn1 * z_regular_over_kr;
            const complex128 No_t_j = complex128(sin_phi, 0.0) * tau_n * dpsi_regular_over_kr;
            const complex128 No_p_j = complex128(cos_phi, 0.0) * pi_n  * dpsi_regular_over_kr;

            if (!is_inside) {

                // Outgoing (scattered) uses outside argument
                const complex128 Mo_r_h(0.0, 0.0);
                const complex128 Mo_t_h = complex128(cos_phi, 0.0) * pi_n  * zn_hankel;
                const complex128 Mo_p_h = complex128(-sin_phi, 0.0) * tau_n * zn_hankel;

                const complex128 Me_r_h(0.0, 0.0);
                const complex128 Me_t_h = complex128(sin_phi, 0.0) * pi_n  * zn_hankel;
                const complex128 Me_p_h = complex128(cos_phi, 0.0) * tau_n * zn_hankel;

                const complex128 Ne_r_h = complex128(cos_phi, 0.0) * nn1 * Pn1 * z_hankel_over_kr;
                const complex128 Ne_t_h = complex128(cos_phi, 0.0) * tau_n * dpsi_hankel_over_kr;
                const complex128 Ne_p_h = complex128(-sin_phi, 0.0) * pi_n  * dpsi_hankel_over_kr;

                const complex128 No_r_h = complex128(sin_phi, 0.0) * nn1 * Pn1 * z_hankel_over_kr;
                const complex128 No_t_h = complex128(sin_phi, 0.0) * tau_n * dpsi_hankel_over_kr;
                const complex128 No_p_h = complex128(cos_phi, 0.0) * pi_n  * dpsi_hankel_over_kr;

                // Total outside = incident + scattered

                // Incident x
                Er_x     += En * ( Mo_r_j - i_unit * Ne_r_j );
                Etheta_x += En * ( Mo_t_j - i_unit * Ne_t_j );
                Ephi_x   += En * ( Mo_p_j - i_unit * Ne_p_j );

                // Scattered x
                Er_x     += En * ( i_unit * a_n * Ne_r_h - b_n * Mo_r_h );
                Etheta_x += En * ( i_unit * a_n * Ne_t_h - b_n * Mo_t_h );
                Ephi_x   += En * ( i_unit * a_n * Ne_p_h - b_n * Mo_p_h );

                // Incident y
                Er_y     += En * ( Me_r_j - i_unit * No_r_j );
                Etheta_y += En * ( Me_t_j - i_unit * No_t_j );
                Ephi_y   += En * ( Me_p_j - i_unit * No_p_j );

                // Scattered y
                Er_y     += En * ( i_unit * a_n * No_r_h - b_n * Me_r_h );
                Etheta_y += En * ( i_unit * a_n * No_t_h - b_n * Me_t_h );
                Ephi_y   += En * ( i_unit * a_n * No_p_h - b_n * Me_p_h );

            } else {

                // Total inside = internal field only
                // Internal field must be evaluated with k_particle r (handled by kr_used_for_regular)

                // Internal x
                Er_x     += En * ( c_n * Mo_r_j - i_unit * d_n * Ne_r_j );
                Etheta_x += En * ( c_n * Mo_t_j - i_unit * d_n * Ne_t_j );
                Ephi_x   += En * ( c_n * Mo_p_j - i_unit * d_n * Ne_p_j );

                // Internal y
                Er_y     += En * ( c_n * Me_r_j - i_unit * d_n * No_r_j );
                Etheta_y += En * ( c_n * Me_t_j - i_unit * d_n * No_t_j );
                Ephi_y   += En * ( c_n * Me_p_j - i_unit * d_n * No_p_j );
            }
        }

        const complex128 Er     = E0x * Er_x     + E0y * Er_y;
        const complex128 Etheta = E0x * Etheta_x + E0y * Etheta_y;
        const complex128 Ephi   = E0x * Ephi_x   + E0y * Ephi_y;

        const complex128 Ex =
            Er * (sin_theta * cos_phi)
            + Etheta * (cos_theta * cos_phi)
            + Ephi * (-sin_phi);

        const complex128 Ey =
            Er * (sin_theta * sin_phi)
            + Etheta * (cos_theta * sin_phi)
            + Ephi * (cos_phi);

        const complex128 Ez =
            Er * (cos_theta)
            + Etheta * (-sin_theta);

        if (field_type == "Ex") {
            field_values[point_index] = Ex;
        } else if (field_type == "Ey") {
            field_values[point_index] = Ey;
        } else if (field_type == "Ez") {
            field_values[point_index] = Ez;
        } else {
            const double mag = std::sqrt(std::norm(Ex) + std::norm(Ey) + std::norm(Ez));
            field_values[point_index] = complex128(mag, 0.0);
        }
    }

    return field_values;
}


std::vector<complex128>
Sphere::compute_scattered_nearfields(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::string& field_type
)
{
    if (field_type != "Ex" && field_type != "Ey" && field_type != "Ez" && field_type != "|E|") {
        throw std::invalid_argument("Invalid field_type. Must be one of: Ex, Ey, Ez, |E|");
    }

    const std::size_t number_of_points = x.size();
    if (number_of_points != y.size() || number_of_points != z.size()) {
        throw std::invalid_argument("x, y, z vectors must have the same length");
    }

    std::vector<complex128> scattered_values(number_of_points);

    const double radius_particle = 0.5 * this->diameter;

    const double medium_refractive_index = this->medium_refractive_index;

    const double k0 = this->source->wavenumber_vacuum;
    const double k_medium = k0 * medium_refractive_index;

    const std::array<complex128, 2> jones_vector = this->source->get_jones_vector_first_row();
    const complex128 E0x = jones_vector[0] * this->source->amplitude;
    const complex128 E0y = jones_vector[1] * this->source->amplitude;

    // Total field from your existing routine
    const std::vector<complex128> total_values = this->compute_nearfields(x, y, z, field_type);

    if (field_type != "|E|") {

        for (std::size_t i = 0; i < number_of_points; ++i) {

            const double x_position = x[i];
            const double y_position = y[i];
            const double z_position = z[i];

            const double r_cartesian = std::sqrt(
                x_position * x_position + y_position * y_position + z_position * z_position
            );

            const bool is_inside = (r_cartesian < radius_particle);

            if (is_inside) {
                // By definition, "scattered" field is typically only defined outside.
                // Keep it explicit.
                scattered_values[i] = complex128(0.0, 0.0);
                continue;
            }

            // Incident plane wave assumed to propagate along +z, consistent with your BH expansion
            const complex128 phase = std::exp(complex128(0.0, 1.0) * complex128(k_medium * z_position, 0.0));

            const complex128 Ex_inc = E0x * phase;
            const complex128 Ey_inc = E0y * phase;
            const complex128 Ez_inc = complex128(0.0, 0.0);

            if (field_type == "Ex") {
                scattered_values[i] = total_values[i] - Ex_inc;
            } else if (field_type == "Ey") {
                scattered_values[i] = total_values[i] - Ey_inc;
            } else {
                scattered_values[i] = total_values[i] - Ez_inc; // Ez_inc = 0
            }
        }

        return scattered_values;
    }

    // Correct scattered magnitude: |E_total - E_inc|
    const std::vector<complex128> total_x = this->compute_nearfields(x, y, z, "Ex");
    const std::vector<complex128> total_y = this->compute_nearfields(x, y, z, "Ey");
    const std::vector<complex128> total_z = this->compute_nearfields(x, y, z, "Ez");

    for (std::size_t i = 0; i < number_of_points; ++i) {

        const double x_position = x[i];
        const double y_position = y[i];
        const double z_position = z[i];

        const double r_cartesian = std::sqrt(
            x_position * x_position + y_position * y_position + z_position * z_position
        );

        const bool is_inside = (r_cartesian < radius_particle);

        if (is_inside) {
            scattered_values[i] = complex128(0.0, 0.0);
            continue;
        }

        const complex128 phase = std::exp(complex128(0.0, 1.0) * complex128(k_medium * z_position, 0.0));

        const complex128 Ex_inc = E0x * phase;
        const complex128 Ey_inc = E0y * phase;
        const complex128 Ez_inc = complex128(0.0, 0.0);

        const complex128 dx = total_x[i] - Ex_inc;
        const complex128 dy = total_y[i] - Ey_inc;
        const complex128 dz = total_z[i] - Ez_inc;

        const double mag = std::sqrt(std::norm(dx) + std::norm(dy) + std::norm(dz));
        scattered_values[i] = complex128(mag, 0.0);
    }

    return scattered_values;
}
