#include "./coreshell.h"


// ---------------------- Constructors ---------------------------------------
CoreShell::CoreShell(
    double _core_diameter,
    double _shell_thickness,
    complex128 _core_refractive_index,
    complex128 _shell_refractive_index,
    double _medium_refractive_index,
    std::shared_ptr<BaseSource> _source,
    size_t _max_order)
:   BaseScatterer(_max_order, std::move(_source), _medium_refractive_index),
    core_diameter(_core_diameter),
    shell_thickness(_shell_thickness),
    core_refractive_index(_core_refractive_index),
    shell_refractive_index(_shell_refractive_index)
{
    this->total_radius = this->core_diameter / 2 + this->shell_thickness;
    this->compute_cross_section();
    this->compute_size_parameter();
    this->max_order = (_max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : _max_order;
    this->apply_medium();
    this->compute_an_bn(this->max_order);
}



// ---------------------- Methods ---------------------------------------

void CoreShell::compute_size_parameter() {
    this->size_parameter = source->wavenumber * this->total_radius * this->medium_refractive_index;
    this->size_parameter_squared = pow(this->size_parameter, 2);
    this->x_shell = source->wavenumber * this->total_radius * this->medium_refractive_index;
    this->x_core = source->wavenumber * this->medium_refractive_index * (0.5 * this->core_diameter);
}
// void CoreShell::compute_size_parameter() {
//     this->size_parameter = source->wavenumber * this->total_radius / 2;// * this->medium_refractive_index;
//     this->size_parameter_squared = pow(this->size_parameter, 2);
//     this->x_shell = source->wavenumber * this->total_radius / 2.0;// * this->medium_refractive_index;
//     this->x_core = source->wavenumber * this->core_diameter / 2.0;// * this->medium_refractive_index;
// }

void CoreShell::compute_cross_section() {
    this->cross_section = Constants::PI * pow(this->total_radius, 2);
}

void CoreShell::apply_medium() {
    this->core_refractive_index /= this->medium_refractive_index;
    this->shell_refractive_index /= this->medium_refractive_index;
    // this->core_diameter *= this->medium_refractive_index;
    // this->shell_thickness *= this->medium_refractive_index;
    // this->total_radius *= this->medium_refractive_index;
}

void CoreShell::compute_an_bn(const size_t max_order)
{
    an.resize(max_order);
    bn.resize(max_order);

    // Calculate scaled parameters and initialize phase shift factors
    complex128
        relative_index = this->shell_refractive_index / this->core_refractive_index,
        u = this->core_refractive_index * this->x_core,
        v = this->shell_refractive_index * this->x_core,
        w = this->shell_refractive_index * this->x_shell,
        sv = sqrt(0.5 * Constants::PI * v),
        sw = sqrt(0.5 * Constants::PI * w),
        sy = sqrt(0.5 * Constants::PI * this->x_shell);

    // Determine the necessary array size for continuity factors
    size_t mx = static_cast<size_t>(std::max( abs( this->core_refractive_index * this->x_shell ), abs( this->shell_refractive_index*this->x_shell ) ));
    size_t nmx  = std::max( max_order, mx ) + 16  ;

    std::vector<complex128> pv(max_order + 1), pw(max_order + 1), py(max_order + 1), chv(max_order + 1), chw(max_order + 1), chy(max_order + 1), gsy(max_order + 1), gs1y(max_order + 1);
    std::vector<complex128> p1y(max_order + 2), ch1y(max_order + 2);

    p1y[0] = sin(this->x_shell);
    ch1y[0] = cos(this->x_shell);

    for (size_t order = 0; order < max_order + 1; order++){
        double nu = order + 1.5 ;

        pw[order] = sw * Cylindrical_::Jn(nu, (complex128) w);
        pv[order] = sv * Cylindrical_::Jn(nu, (complex128) v);
        py[order] = sy * Cylindrical_::Jn(nu, (complex128) this->x_shell);

        chv[order] = -sv * Cylindrical_::Yn(nu, v);
        chw[order] = -sw * Cylindrical_::Yn(nu, w);
        chy[order] = -sy * Cylindrical_::Yn(nu, this->x_shell);

        p1y[order + 1] = py[order];
        ch1y[order + 1] = chy[order];
        gsy[order] = py[order]  - complex128(0, 1) * chy[order];
        gs1y[order] = p1y[order] - complex128(0, 1) * ch1y[order];
    }

    // Calculate continuity factors in reverse order
    std::vector<complex128> Du(nmx, 0.0), Dv(nmx, 0.0), Dw(nmx, 0.0);

    for (int i = nmx - 1; i > 1; i--){
        Du[i-1] = (double)i / u - 1.0 / (Du[i] + (double)i / u);
        Dv[i-1] = (double)i / v - 1.0 / (Dv[i] + (double)i / v);
        Dw[i-1] = (double)i / w - 1.0 / (Dw[i] + (double)i / w);
    }

    // Resize continuity factors to maximum order needed
    Du.erase(Du.begin());
    Dv.erase(Dv.begin());
    Dw.erase(Dw.begin());

    // Calculate Mie coefficients
    std::vector<complex128> uu(max_order), vv(max_order), fv(max_order), dns(max_order), gns(max_order), a1(max_order), b1(max_order);

    for (size_t order=0; order < max_order; order++){
        double idx = static_cast<double>(order + 1);

        uu[order] = relative_index * Du[order] - Dv[order];
        vv[order] = Du[order] / relative_index - Dv[order];
        fv[order] = pv[order] / chv[order];
        dns[order] = ((uu[order] * fv[order] / pw[order]) / (uu[order] * (pw[order] - chw[order] * fv[order]) + pw[order] / pv[order] / chv[order])) + Dw[order];
        gns[order] = ((vv[order] * fv[order] / pw[order]) / (vv[order] * (pw[order] - chw[order] * fv[order]) + pw[order] / pv[order] / chv[order])) + Dw[order];
        a1[order] = dns[order] / shell_refractive_index + idx / x_shell;
        b1[order] = shell_refractive_index * gns[order] + idx / x_shell;
        an[order] = (py[order] * a1[order] - p1y[order]) / (gsy[order] * a1[order] - gs1y[order]);
        bn[order] = (py[order] * b1[order] - p1y[order]) / (gsy[order] * b1[order] - gs1y[order]);
    }
}


double CoreShell::get_Qsca() const {
    double value = 0;

    for(size_t it = 0; it < max_order; ++it){
        double n = (double) it + 1;
        value += (2. * n + 1.) * ( pow( std::abs(this->an[it]), 2) + pow( std::abs(this->bn[it]), 2)  );
    }
    return value * 2. / size_parameter_squared;
}

double CoreShell::get_Qext() const {
    double value = 0;
    for(size_t it = 0; it < max_order; ++it)
    {
        double n = (double) it + 1;
        value += (2.* n + 1.) * std::real( this->an[it] + this->bn[it] );

    }
    return value * 2. / size_parameter_squared;
}

double CoreShell::get_Qback() const {
    complex128 value = 0;

    for(size_t it = 0; it < max_order-1; ++it)
    {
        double n = (double) it + 1;

        value += (2. * n + 1) * pow(-1., n) * ( this->an[it] - this->bn[it] ) ;
    }

    value = pow( std::abs(value), 2. ) / size_parameter_squared;
    return std::abs(value);
}

double CoreShell::get_Qforward() const {
    complex128 value = 0;

    for(size_t it = 0; it < max_order-1; ++it)
    {
        double n = (double) it + 1;

        value += (2. * n + 1) * ( this->an[it] + this->bn[it] ) ;
    }

    value = pow( std::abs(value), 2. ) / size_parameter_squared;
    return std::abs(value);
}

double CoreShell::get_g() const {
    double value = 0;

      for(size_t it = 0; it < max_order-1; ++it) {
         double n = (double) it + 1;

          value += ( n * (n + 2.) / (n + 1.) ) * std::real(this->an[it] * std::conj(this->an[it+1]) + this->bn[it] * std::conj(this->bn[it+1]) );
          value += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( this->an[it] * std::conj(this->bn[it]) );
      }
      return value * 4. / ( get_Qsca() * size_parameter_squared );
}

std::tuple<std::vector<complex128>, std::vector<complex128>>
CoreShell::compute_s1s2(const std::vector<double> &phi) const {
    std::vector<complex128> S1, S2;

    S1.reserve(phi.size());
    S2.reserve(phi.size());

    std::vector<double> mu, prefactor = get_prefactor();

    mu.reserve(phi.size());

    for (const double phi : phi)
        mu.push_back( cos( phi - Constants::PI / 2.0 ) );

    for (size_t i = 0; i < phi.size(); i++){
        auto [pin, taun] = this->get_pi_tau(mu[i], max_order);
        complex128 S1_temp = 0., S2_temp = 0.;

        for (size_t m = 0; m < max_order ; m++){
            S1_temp += prefactor[m] * ( this->an[m] * pin[m] +  this->bn[m] * taun[m] );
            S2_temp += prefactor[m] * ( this->an[m] * taun[m] + this->bn[m] * pin[m]  );
        }
        S1.push_back(S1_temp);
        S2.push_back(S2_temp);
    }

    return std::make_tuple(std::move(S1), std::move(S2));
}

void CoreShell::compute_cn_dn(size_t) {
    throw std::runtime_error("No implementation of cn and dn exists yet!");
}

std::vector<complex128> CoreShell::compute_nearfields(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::string&) {
    throw std::runtime_error("Near-field computation is not implemented for CoreShell scatterers.");
}

// -
