#include "./coreshell.h"

// ---------------------- Methods ---------------------------------------

void CoreShell::init(const std::shared_ptr<BaseSource>& source, const size_t _max_order) {
    this->core_material->initialize(source->wavelength);
    this->shell_material->initialize(source->wavelength);
    this->medium->initialize(source->wavelength);

    this->total_diameter = this->core_diameter + 2.0 * this->shell_thickness;

    this->compute_cross_section();
    this->compute_size_parameter(source);

    this->max_order = (_max_order == 0)
        ? this->get_wiscombe_criterion(this->size_parameter)
        : _max_order;

    this->compute_an_bn(this->max_order);
}


void CoreShell::compute_size_parameter(const std::shared_ptr<BaseSource>& source) {
    const double medium_refractive_index = this->medium->get_refractive_index();

    this->x_shell = source->wavenumber_vacuum * this->total_diameter / 2.0 * medium_refractive_index;
    this->x_core = source->wavenumber_vacuum * this->core_diameter / 2.0 * medium_refractive_index;

    this->size_parameter = this->x_shell;
    this->size_parameter_squared = std::pow(this->size_parameter, 2);
}

void CoreShell::compute_cross_section() {
    this->cross_section = Constants::PI * std::pow(this->total_diameter / 2.0, 2);
}

void CoreShell::apply_medium() {
    throw std::runtime_error(
        "CoreShell::apply_medium should not be used. "
        "Relative refractive indices are computed locally in compute_an_bn."
    );
}

void CoreShell::compute_an_bn(const size_t max_order)
{
    this->an.resize(max_order);
    this->bn.resize(max_order);

    const complex128 core_refractive_index = this->core_material->get_refractive_index();
    const complex128 shell_refractive_index = this->shell_material->get_refractive_index();
    const double medium_refractive_index = this->medium->get_refractive_index();

    const complex128 core_relative_refractive_index = core_refractive_index / medium_refractive_index;
    const complex128 shell_relative_refractive_index = shell_refractive_index / medium_refractive_index;
    const complex128 shell_to_core_relative_index = shell_relative_refractive_index / core_relative_refractive_index;

    const complex128 u = core_relative_refractive_index * this->x_core;
    const complex128 v = shell_relative_refractive_index * this->x_core;
    const complex128 w = shell_relative_refractive_index * this->x_shell;

    const complex128 sqrt_half_pi_v = std::sqrt(0.5 * Constants::PI * v);
    const complex128 sqrt_half_pi_w = std::sqrt(0.5 * Constants::PI * w);
    const complex128 sqrt_half_pi_x_shell = std::sqrt(0.5 * Constants::PI * this->x_shell);

    const size_t max_argument_index = static_cast<size_t>(
        std::max(
            std::abs(core_relative_refractive_index * this->x_shell),
            std::abs(shell_relative_refractive_index * this->x_shell)
        )
    );

    const size_t number_of_downward_terms = std::max(max_order, max_argument_index) + 16;

    std::vector<complex128> pv(max_order + 1);
    std::vector<complex128> pw(max_order + 1);
    std::vector<complex128> py(max_order + 1);

    std::vector<complex128> chv(max_order + 1);
    std::vector<complex128> chw(max_order + 1);
    std::vector<complex128> chy(max_order + 1);

    std::vector<complex128> gsy(max_order + 1);
    std::vector<complex128> gs1y(max_order + 1);

    std::vector<complex128> p1y(max_order + 2);
    std::vector<complex128> ch1y(max_order + 2);

    p1y[0] = std::sin(this->x_shell);
    ch1y[0] = std::cos(this->x_shell);

    for (size_t order = 0; order < max_order + 1; ++order) {
        const double nu = static_cast<double>(order) + 1.5;

        pw[order] = sqrt_half_pi_w * Cylindrical_::Jn(nu, w);
        pv[order] = sqrt_half_pi_v * Cylindrical_::Jn(nu, v);
        py[order] = sqrt_half_pi_x_shell * Cylindrical_::Jn(nu, this->x_shell);

        chv[order] = -sqrt_half_pi_v * Cylindrical_::Yn(nu, v);
        chw[order] = -sqrt_half_pi_w * Cylindrical_::Yn(nu, w);
        chy[order] = -sqrt_half_pi_x_shell * Cylindrical_::Yn(nu, this->x_shell);

        p1y[order + 1] = py[order];
        ch1y[order + 1] = chy[order];

        gsy[order] = py[order] - complex128(0.0, 1.0) * chy[order];
        gs1y[order] = p1y[order] - complex128(0.0, 1.0) * ch1y[order];
    }

    std::vector<complex128> Du(number_of_downward_terms, 0.0);
    std::vector<complex128> Dv(number_of_downward_terms, 0.0);
    std::vector<complex128> Dw(number_of_downward_terms, 0.0);

    for (int index = static_cast<int>(number_of_downward_terms) - 1; index > 1; --index) {
        Du[index - 1] = static_cast<double>(index) / u - 1.0 / (Du[index] + static_cast<double>(index) / u);
        Dv[index - 1] = static_cast<double>(index) / v - 1.0 / (Dv[index] + static_cast<double>(index) / v);
        Dw[index - 1] = static_cast<double>(index) / w - 1.0 / (Dw[index] + static_cast<double>(index) / w);
    }

    Du.erase(Du.begin());
    Dv.erase(Dv.begin());
    Dw.erase(Dw.begin());

    std::vector<complex128> uu(max_order);
    std::vector<complex128> vv(max_order);
    std::vector<complex128> fv(max_order);
    std::vector<complex128> dns(max_order);
    std::vector<complex128> gns(max_order);
    std::vector<complex128> a1(max_order);
    std::vector<complex128> b1(max_order);

    for (size_t order = 0; order < max_order; ++order) {
        const double idx = static_cast<double>(order + 1);

        uu[order] = shell_to_core_relative_index * Du[order] - Dv[order];
        vv[order] = Du[order] / shell_to_core_relative_index - Dv[order];
        fv[order] = pv[order] / chv[order];

        dns[order] =
            (
                (uu[order] * fv[order] / pw[order]) /
                (
                    uu[order] * (pw[order] - chw[order] * fv[order]) +
                    pw[order] / pv[order] / chv[order]
                )
            ) + Dw[order];

        gns[order] =
            (
                (vv[order] * fv[order] / pw[order]) /
                (
                    vv[order] * (pw[order] - chw[order] * fv[order]) +
                    pw[order] / pv[order] / chv[order]
                )
            ) + Dw[order];

        a1[order] = dns[order] / shell_relative_refractive_index + idx / this->x_shell;
        b1[order] = shell_relative_refractive_index * gns[order] + idx / this->x_shell;

        this->an[order] = (py[order] * a1[order] - p1y[order]) / (gsy[order] * a1[order] - gs1y[order]);
        this->bn[order] = (py[order] * b1[order] - p1y[order]) / (gsy[order] * b1[order] - gs1y[order]);
    }
}


double CoreShell::get_Qsca() const {
    double value = 0.0;

    for (size_t it = 0; it < this->max_order; ++it) {
        const double n = static_cast<double>(it) + 1.0;
        value += (2.0 * n + 1.0) * (std::pow(std::abs(this->an[it]), 2) + std::pow(std::abs(this->bn[it]), 2));
    }

    return value * 2.0 / this->size_parameter_squared;
}

double CoreShell::get_Qext() const {
    double value = 0.0;

    for (size_t it = 0; it < this->max_order; ++it) {
        const double n = static_cast<double>(it) + 1.0;
        value += (2.0 * n + 1.0) * std::real(this->an[it] + this->bn[it]);
    }

    return value * 2.0 / this->size_parameter_squared;
}

double CoreShell::get_Qback() const {
    complex128 value = 0.0;

    for (size_t it = 0; it < this->max_order; ++it) {
        const double n = static_cast<double>(it) + 1.0;
        value += (2.0 * n + 1.0) * std::pow(-1.0, n) * (this->an[it] - this->bn[it]);
    }

    value = std::pow(std::abs(value), 2.0) / this->size_parameter_squared;
    return std::abs(value);
}

double CoreShell::get_Qforward() const {
    complex128 value = 0.0;

    for (size_t it = 0; it < this->max_order; ++it) {
        const double n = static_cast<double>(it) + 1.0;
        value += (2.0 * n + 1.0) * (this->an[it] + this->bn[it]);
    }

    value = std::pow(std::abs(value), 2.0) / this->size_parameter_squared;
    return std::abs(value);
}

double CoreShell::get_g() const {
    double value = 0.0;

    for (size_t it = 0; it < this->max_order - 1; ++it) {
        const double n = static_cast<double>(it) + 1.0;

        value +=
            (n * (n + 2.0) / (n + 1.0)) *
            std::real(
                this->an[it] * std::conj(this->an[it + 1]) +
                this->bn[it] * std::conj(this->bn[it + 1])
            );

        value +=
            ((2.0 * n + 1.0) / (n * (n + 1.0))) *
            std::real(this->an[it] * std::conj(this->bn[it]));
    }

    return value * 4.0 / (this->get_Qsca() * this->size_parameter_squared);
}

std::pair<std::vector<complex128>, std::vector<complex128>>
CoreShell::compute_s1s2(const std::vector<double>& phi) const {
    std::vector<complex128> S1;
    std::vector<complex128> S2;

    S1.reserve(phi.size());
    S2.reserve(phi.size());

    std::vector<double> mu;
    mu.reserve(phi.size());

    const std::vector<double> prefactor = this->get_prefactor();

    for (const double angle : phi) {
        mu.push_back(std::cos(angle - Constants::PI / 2.0));
    }

    for (size_t index = 0; index < phi.size(); ++index) {
        auto [pin, taun] = this->get_pi_tau(mu[index], this->max_order);

        complex128 S1_temp = 0.0;
        complex128 S2_temp = 0.0;

        for (size_t order = 0; order < this->max_order; ++order) {
            S1_temp += prefactor[order] * (this->an[order] * pin[order] + this->bn[order] * taun[order]);
            S2_temp += prefactor[order] * (this->an[order] * taun[order] + this->bn[order] * pin[order]);
        }

        S1.push_back(S1_temp);
        S2.push_back(S2_temp);
    }

    return std::make_pair(std::move(S1), std::move(S2));
}

void CoreShell::compute_cn_dn(size_t) {
    throw std::runtime_error("No implementation of cn and dn exists yet!");
}