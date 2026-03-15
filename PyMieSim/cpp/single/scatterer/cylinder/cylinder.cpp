#include "./cylinder.h"


// ---------------------- Methods ---------------------------------------

void InfiniteCylinder::init(const std::shared_ptr<BaseSource>& source, const size_t _max_order) {
    this->material->initialize(source->wavelength);
    this->medium->initialize(source->wavelength);
    this->jones_vector = source->get_jones_vector_first_row();

    this->compute_cross_section();
    this->compute_size_parameter(source);
    this->max_order = (_max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : _max_order;
    this->compute_an_bn(this->max_order);
}


void InfiniteCylinder::compute_size_parameter(const std::shared_ptr<BaseSource>& source) {
    this->size_parameter = source->wavenumber_vacuum * this->diameter / 2;
    this->size_parameter_squared = pow(this->size_parameter, 2);
}

void InfiniteCylinder::compute_cross_section() {
    this->cross_section = this->diameter;
}

double InfiniteCylinder::get_Qsca() const {
    double Qsca1 = 0.0;
    double Qsca2 = 0.0;

    for (size_t order = 1; order < max_order; ++order) {
        Qsca1 += std::norm(this->a1n[order]) + std::norm(this->b1n[order]);
        Qsca2 += std::norm(this->a2n[order]) + std::norm(this->b2n[order]);
    }

    Qsca1 = 2.0 / size_parameter * (2.0 * Qsca1 + std::norm(this->b1n[0]));
    Qsca2 = 2.0 / size_parameter * (2.0 * Qsca2 + std::norm(this->a2n[0]));

    return this->process_polarization(Qsca1, Qsca2);
}

double InfiniteCylinder::get_Qext() const {
    double Qext1_sum = 0.0;
    double Qext2_sum = 0.0;

    for (size_t order = 1; order < max_order; ++order) {
        Qext1_sum += std::real(this->b1n[order]);
        Qext2_sum += std::real(this->a2n[order]);
    }

    const double Qext1 = 2.0 / size_parameter * (std::real(this->b1n[0]) + 2.0 * Qext1_sum);
    const double Qext2 = 2.0 / size_parameter * (std::real(this->a2n[0]) + 2.0 * Qext2_sum);

    return this->process_polarization(Qext1, Qext2);
}

double InfiniteCylinder::process_polarization(const double value_0, const double value_1) const {
    const double Ex_weight = std::norm(this->jones_vector[0]);
    const double Ey_weight = std::norm(this->jones_vector[1]);

    return value_1 * Ex_weight + value_0 * Ey_weight;
}

void InfiniteCylinder::compute_an_bn(const size_t max_order) {
    // Resize vectors to hold Mie coefficients for the specified maximum order
    this->a1n.resize(max_order);
    this->b1n.resize(max_order);
    this->a2n.resize(max_order);
    this->b2n.resize(max_order);

    // Compute the arguments used in the Bessel and Hankel function calculations
    complex128
        m = this->material->get_refractive_index() / this->medium->get_refractive_index(), // Relative refractive index
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

void InfiniteCylinder::compute_cn_dn(const size_t max_order) {
    // Resize vectors to hold internal field coefficients for the specified maximum order
    this->c1n.resize(max_order);
    this->d1n.resize(max_order);
    this->c2n.resize(max_order);
    this->d2n.resize(max_order);

    // Compute the relative refractive index and scaled internal size parameter
    complex128 m = this->material->get_refractive_index() / this->medium->get_refractive_index();
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


std::pair<std::vector<complex128>, std::vector<complex128>>
InfiniteCylinder::compute_s1s2(const std::vector<double> &phi) const{
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
InfiniteCylinder::compute_dn(double nmx, complex128 z) const { //Page 205 of BH
    std::vector<complex128> Dn(nmx, 0.0);

    for (double n = nmx - 1; n > 0; n--)
        Dn[n-1] = n/z - ( 1. / (Dn[n] + n/z) );

    return Dn;
}
