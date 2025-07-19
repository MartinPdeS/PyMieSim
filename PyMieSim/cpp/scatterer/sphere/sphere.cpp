#include "scatterer/sphere/sphere.h"


// ---------------------- Constructors ---------------------------------------
Sphere::Sphere(const double diameter, const complex128 refractive_index, const double medium_refractive_index, const BaseSource &source, size_t max_order, bool compute_cn_dn)
: BaseScatterer(max_order, source, medium_refractive_index), diameter(diameter), refractive_index(refractive_index)
{
    this->compute_cross_section();
    this->compute_size_parameter();
    this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
    this->compute_an_bn(this->max_order);

    if (compute_cn_dn)
        this->compute_cn_dn(this->max_order);
}

// ---------------------- Methods ---------------------------------------
void Sphere::compute_size_parameter() {
    this->size_parameter = source.wavenumber * this->diameter / 2 * this->medium_refractive_index;
    this->size_parameter_squared = pow(this->size_parameter, 2);
}

void Sphere::compute_cross_section() {
    this->cross_section = PI * std::pow(this->diameter / 2.0, 2);
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
        psi_n = +size_parameter * compute_jn(nu, size_parameter);
        chi_n = -size_parameter * compute_yn(nu, size_parameter);

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
        m = this->refractive_index,
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
        jnx.push_back(compute_jn(order + 1, x));

        jnmx.push_back(1.0 / (compute_jn(order + 1, z)));
        yx.push_back(compute_yn(order + 1, x));
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
        mu.push_back( cos( phi - PI / 2.0 ) );

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

// -
