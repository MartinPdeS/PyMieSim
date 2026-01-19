#include "./cylinder.h"

// ---------------------- Constructors ---------------------------------------
Cylinder::Cylinder(double _diameter, complex128 _refractive_index, double _medium_refractive_index, std::shared_ptr<BaseSource> _source, size_t _max_order) :
BaseScatterer(_max_order, std::move(_source), _medium_refractive_index), diameter(_diameter), refractive_index(_refractive_index)
{
    this->compute_cross_section();
    this->compute_size_parameter();
    this->max_order = (_max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : _max_order;
    this->compute_an_bn(this->max_order);
}

// ---------------------- Methods ---------------------------------------

void Cylinder::compute_size_parameter() {
    this->size_parameter = source->wavenumber_vacuum * this->diameter / 2;
    this->size_parameter_squared = pow(this->size_parameter, 2);
}

void Cylinder::compute_cross_section() {
    this->cross_section = this->diameter;
}

double Cylinder::get_Qsca() const {
    complex128 Qsca1=0, Qsca2=0;

    for(size_t order = 1; order < max_order; order++)
    {
        Qsca1 +=  pow( std::abs(this->a1n[order]), 2 ) + pow( std::abs(this->b1n[order]), 2 ) ;
        Qsca2 +=  pow( std::abs(this->a2n[order]), 2 ) + pow( std::abs(this->b2n[order]), 2 ) ;
    }

    Qsca1 =  2. / size_parameter * ( 2.0 * Qsca1 + pow( std::abs(this->b1n[0]), 2 ) );
    Qsca2 =  2. / size_parameter * ( 2.0 * Qsca2 + pow( std::abs(this->a2n[0]), 2 ) );

    return process_polarization(Qsca1, Qsca2);
}

double Cylinder::get_Qext() const {
    complex128 Qext1 = 0, Qext2 = 0;

    for(size_t it = 1; it < max_order; ++it){
        Qext1 += this->b1n[it];
        Qext2 += this->a2n[it];
    }

    Qext1 = 2. / size_parameter * std::real( this->b1n[0] + 2.0 * Qext1 );
    Qext2 = 2. / size_parameter * std::real( this->a1n[0] + 2.0 * Qext2 );

    return this->process_polarization(Qext1, Qext2);
}

double Cylinder::get_g() const {
    return this->get_g_with_farfields(1000);
}

double Cylinder::process_polarization(const complex128 value_0, const complex128 value_1) const {
    std::array<complex128, 2> jones_vector = this->source->get_jones_vector_first_row();
    complex128
        Ex = jones_vector[0],
        Ey = jones_vector[1];

    return std::abs( value_1 ) * pow(std::abs(Ex), 2) + std::abs( value_0 ) * pow(std::abs(Ey), 2);
}

void Cylinder::compute_an_bn(const size_t max_order) {
    // Resize vectors to hold Mie coefficients for the specified maximum order
    this->a1n.resize(max_order);
    this->b1n.resize(max_order);
    this->a2n.resize(max_order);
    this->b2n.resize(max_order);

    // Compute the arguments used in the Bessel and Hankel function calculations
    complex128
        m = this->refractive_index / this->medium_refractive_index, // Relative refractive index
        mx = m * size_parameter; // Scaled size parameter for internal calculations

    // Precompute Bessel and Hankel functions up to max_order
    std::vector<complex128>
        J_z(max_order + 1),
        J_z_p(max_order + 1),
        J_x(max_order + 1),
        J_x_p(max_order + 1),
        H_x(max_order + 1),
        H_x_p(max_order + 1);

    for (size_t order = 0; order < max_order + 1; ++order){
        J_z[order] = Cylindrical_::Jn(order, mx);
        J_z_p[order] = Cylindrical_::Jnp(order, mx);
        J_x[order] = Cylindrical_::Jn(order, size_parameter);
        J_x_p[order] = Cylindrical_::Jnp(order, size_parameter);
        H_x[order] = Cylindrical_::H1n(order, size_parameter);
        H_x_p[order] = Cylindrical_::H1np(order, size_parameter);
    }

    // Compute Mie coefficients a1n, a2n, b1n, b2n for each order
    for (size_t order = 0; order < max_order; order++){
        complex128 numerator_a = m * J_z[order] * J_x_p[order] - J_z_p[order] * J_x[order];
        complex128 denominator_a = m * J_z[order] * H_x_p[order] - J_z_p[order] * H_x[order];
        this->a1n[order] = 0.0 ;
        this->a2n[order] = numerator_a / denominator_a;

        complex128 numerator_b = J_z[order] * J_x_p[order] - m * J_z_p[order] * J_x[order];
        complex128 denominator_b = J_z[order] * H_x_p[order] - m * J_z_p[order] * H_x[order];
        this->b1n[order] = numerator_b / denominator_b;
        this->b2n[order] = 0.0 ;
    }
}

void Cylinder::compute_cn_dn(const size_t max_order) {
    // Resize vectors to hold internal field coefficients for the specified maximum order
    this->c1n.resize(max_order);
    this->d1n.resize(max_order);
    this->c2n.resize(max_order);
    this->d2n.resize(max_order);

    // Compute the relative refractive index and scaled internal size parameter
    complex128 m = this->refractive_index / this->medium_refractive_index;
    complex128 mx = m * size_parameter;

    // Impedance ratio (used in TE mode internal coefficient)
    // complex128 m_tilde = this->medium_impedance / this->particle_impedance;

    // Precompute Bessel and Hankel functions
    std::vector<complex128>
        J_z(max_order + 1),
        J_z_p(max_order + 1),
        J_x(max_order + 1),
        J_x_p(max_order + 1),
        H_x(max_order + 1),
        H_x_p(max_order + 1);

    for (size_t order = 0; order < max_order + 1; ++order){
        J_z[order]  = Cylindrical_::Jn(order, mx);
        J_z_p[order]= Cylindrical_::Jnp(order, mx);
        J_x[order]  = Cylindrical_::Jn(order, size_parameter);
        J_x_p[order]= Cylindrical_::Jnp(order, size_parameter);
        H_x[order]  = Cylindrical_::H1n(order, size_parameter);
        H_x_p[order]= Cylindrical_::H1np(order, size_parameter);
    }

    // Compute internal coefficients c1n, c2n, d1n, d2n
    for (size_t order = 0; order < max_order; ++order) {
        // --- TM polarization (Mode I) ---
        complex128 num_d1 = J_z[order] * H_x_p[order] - m * J_z_p[order] * H_x[order];
        complex128 den_d1 = J_z[order] * J_x_p[order] - m * J_z_p[order] * J_x[order];
        this->d1n[order] = num_d1 / den_d1;
        this->c1n[order] = 0.0;  // Not used in TM polarization

        // --- TE polarization (Mode II) ---
        complex128 num_c2 = m * J_z[order] * H_x_p[order] - J_z_p[order] * H_x[order];
        complex128 den_c2 = m * J_z[order] * J_x_p[order] - J_z_p[order] * J_x[order];
        this->c2n[order] = num_c2 / den_c2;
        this->d2n[order] = 0.0;  // Not used in TE polarization
    }
}


std::tuple<std::vector<complex128>, std::vector<complex128>>
Cylinder::compute_s1s2(const std::vector<double> &phi) const{
    std::vector<complex128> T1(phi.size()), T2(phi.size());

    for (size_t i = 0; i < phi.size(); i++){
        T1[i] = this->b1n[0];
        T2[i] = this->a2n[0];
        for (size_t order = 1; order < max_order ; order++){
            T1[i] += 2.0 * this->b1n[order] * cos(order * (Constants::PI - (phi[i] + Constants::PI/2.0) ) );
            T2[i] += 2.0 * this->a2n[order] * cos(order * (Constants::PI - (phi[i] + Constants::PI/2.0) ) );
        }
    }

    return std::make_tuple(std::move(T1), std::move(T2));
}


std::vector<complex128>
Cylinder::compute_dn(double nmx, complex128 z) const { //Page 205 of BH
    std::vector<complex128> Dn(nmx, 0.0);

    for (double n = nmx - 1; n > 0; n--)
        Dn[n-1] = n/z - ( 1. / (Dn[n] + n/z) );

    return Dn;
}

std::vector<complex128> Cylinder::compute_nearfields(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::string&) {
    throw std::runtime_error("Near-field computation is not implemented for Cylinder scatterers.");
}

// std::vector<complex128>
// Cylinder::compute_near_field(const std::vector<double>& x,
//                              const std::vector<double>& y,
//                              const std::string& field_type)
// {
//     if (cn.empty() || dn.empty())
//         throw std::runtime_error("Near-field computation requires cn and dn coefficients.");

//     const size_t n_points = x.size();
//     if (n_points != y.size())
//         throw std::invalid_argument("x and y vectors must have the same length");

//     std::vector<complex128> field_values(n_points);
//     const double k = source->wavenumber_vacuum * medium_refractive_index;
//     const complex128 i(0.0, 1.0);

//     for (size_t point_idx = 0; point_idx < n_points; ++point_idx) {
//         const double x_pos = x[point_idx];
//         const double y_pos = y[point_idx];

//         const double r = std::hypot(x_pos, y_pos);
//         const double phi = std::atan2(y_pos, x_pos);
//         const bool is_inside = (r < this->radius);

//         complex128 Ez = 0.0;

//         for (int n = -static_cast<int>(max_order); n <= static_cast<int>(max_order); ++n) {
//             size_t idx = std::abs(n);
//             complex128 coeff = is_inside ? dn[idx] : bn[idx];

//             double arg = k * r;
//             complex128 radial_func;

//             if (is_inside)
//                 radial_func = cyl_bessel_j(n, arg);
//             else
//                 radial_func = cyl_hankel_1(n, arg);

//             complex128 exp_term = std::exp(i * static_cast<double>(n) * phi);
//             Ez += coeff * radial_func * exp_term;
//         }

//         if (field_type == "Ez") {
//             field_values[point_idx] = Ez;
//         } else if (field_type == "|E|") {
//             field_values[point_idx] = std::abs(Ez); // scalar Ez only, no transverse components in this simplified TM case
//         } else {
//             throw std::invalid_argument("Unsupported field_type for cylinder: must be Ez or |E| (TM mode only in current implementation)");
//         }
//     }

//     return field_values;
// }
