#include "./coreshell.h"
#include <array>
#include <functional>
#include <tuple>

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

namespace {
using C = complex128;

struct RadialFunctions {
    C z;
    C y;
    C derivative;
};

RadialFunctions radial_functions(const size_t n, const C& kr) {
    const C j = Spherical_::jn(static_cast<double>(n), kr);
    const C y = Spherical_::yn(static_cast<double>(n), kr);
    const C jm1 = Spherical_::jn(static_cast<double>(n - 1), kr);
    return {j, y, (kr * jm1 - static_cast<double>(n) * j) / kr};
}

RadialFunctions outgoing_functions(const size_t n, const C& kr) {
    const C h = Spherical_::H1n(static_cast<double>(n), kr);
    const C hm1 = Spherical_::H1n(static_cast<double>(n - 1), kr);
    return {h, C{0.0, 0.0}, (kr * hm1 - static_cast<double>(n) * h) / kr};
}

std::array<C, 2> solve_two_by_two(const C& a, const C& b, const C& c, const C& d,
                                   const C& rhs_a, const C& rhs_b) {
    const C determinant = a * d - b * c;
    return {(rhs_a * d - b * rhs_b) / determinant,
            (a * rhs_b - rhs_a * c) / determinant};
}
}

void CoreShell::compute_cn_dn(const size_t requested_order) {
    const size_t order_count = requested_order == 0 ? this->max_order : requested_order;
    this->cn.resize(order_count);
    this->dn.resize(order_count);
    this->shell_m_regular.resize(order_count);
    this->shell_m_irregular.resize(order_count);
    this->shell_n_regular.resize(order_count);
    this->shell_n_irregular.resize(order_count);

    const C medium_index = this->medium->get_refractive_index();
    const C core_index = this->core_material->get_refractive_index() / medium_index;
    const C shell_index = this->shell_material->get_refractive_index() / medium_index;
    const C core_radius_argument = core_index * this->x_core;
    const C shell_inner_argument = shell_index * this->x_core;
    const C shell_outer_argument = shell_index * this->x_shell;
    const C outer_argument = C{this->x_shell, 0.0};

    // Solve the two interface conditions independently for the TE (M) and
    // TM (N) multipoles.  This is a little more verbose than the closed
    // Aden–Kerker expressions, but is much easier to audit and handles the
    // three radial regions explicitly.
    for (size_t index = 0; index < order_count; ++index) {
        const size_t n = index + 1;
        const auto shell_j = radial_functions(n, shell_outer_argument);
        const auto shell_i = radial_functions(n, shell_inner_argument);
        const auto outer_j = radial_functions(n, outer_argument);
        const auto outer_h = outgoing_functions(n, outer_argument);

        const C outer_m = outer_j.z - this->bn[index] * outer_h.z;
        const C outer_m_h = outer_j.derivative - this->bn[index] * outer_h.derivative;
        const auto shell_m = solve_two_by_two(
            shell_j.z, Spherical_::yn(static_cast<double>(n), shell_outer_argument),
            shell_index * shell_j.derivative,
            shell_index * (shell_outer_argument * Spherical_::yn(static_cast<double>(n - 1), shell_outer_argument)
                - static_cast<double>(n) * Spherical_::yn(static_cast<double>(n), shell_outer_argument)) / shell_outer_argument,
            outer_m, outer_m_h
        );

        const C outer_n = outer_j.derivative - this->an[index] * outer_h.derivative;
        const C outer_n_h = outer_j.z - this->an[index] * outer_h.z;
        const auto shell_n = solve_two_by_two(
            shell_j.derivative,
            (shell_outer_argument * Spherical_::yn(static_cast<double>(n - 1), shell_outer_argument)
                - static_cast<double>(n) * Spherical_::yn(static_cast<double>(n), shell_outer_argument)) / shell_outer_argument,
            shell_index * shell_j.z,
            shell_index * Spherical_::yn(static_cast<double>(n), shell_outer_argument),
            outer_n, outer_n_h
        );

        this->shell_m_regular[index] = shell_m[0];
        this->shell_m_irregular[index] = shell_m[1];
        this->shell_n_regular[index] = shell_n[0];
        this->shell_n_irregular[index] = shell_n[1];

        const C core_m_boundary = shell_m[0] * shell_i.z + shell_m[1] * Spherical_::yn(static_cast<double>(n), shell_inner_argument);
        const C core_n_boundary = shell_n[0] * shell_i.derivative + shell_n[1] *
            (shell_inner_argument * Spherical_::yn(static_cast<double>(n - 1), shell_inner_argument)
                - static_cast<double>(n) * Spherical_::yn(static_cast<double>(n), shell_inner_argument)) / shell_inner_argument;

        this->cn[index] = core_m_boundary / Spherical_::jn(static_cast<double>(n), core_radius_argument);
        this->dn[index] = core_n_boundary /
            ((core_radius_argument * Spherical_::jn(static_cast<double>(n - 1), core_radius_argument)
                - static_cast<double>(n) * Spherical_::jn(static_cast<double>(n), core_radius_argument)) / core_radius_argument);
    }
}

namespace {
std::vector<complex128> evaluate_core_shell_field(
    const CoreShell& scatterer,
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::string& field_type,
    const std::shared_ptr<BaseSource>& source,
    const bool scattered_only
) {
    if (field_type != "Ex" && field_type != "Ey" && field_type != "Ez" && field_type != "|E|")
        throw std::invalid_argument("Invalid field_type. Must be one of: Ex, Ey, Ez, |E|");
    if (x.size() != y.size() || x.size() != z.size())
        throw std::invalid_argument("x, y, z vectors must have the same length");

    const C i_unit{0.0, 1.0};
    const C e0x = source->polarization.jones_vector[0] * source->amplitude;
    const C e0y = source->polarization.jones_vector[1] * source->amplitude;
    const double k_medium = source->wavenumber_vacuum * scatterer.medium->get_refractive_index();
    const double core_radius = scatterer.core_diameter / 2.0;
    const double outer_radius = scatterer.total_diameter / 2.0;
    const C core_index = scatterer.core_material->get_refractive_index();
    const C shell_index = scatterer.shell_material->get_refractive_index();
    std::vector<C> values(x.size(), C{0.0, 0.0});

    auto clamp = [](double value) { return std::max(-1.0, std::min(1.0, value)); };
    auto radial = [](size_t n, const C& argument, int kind) {
        if (kind == 2) return outgoing_functions(n, argument);
        auto result = radial_functions(n, argument);
        if (kind == 1) {
            result.z = result.y;
            result.derivative = (argument * Spherical_::yn(static_cast<double>(n - 1), argument)
                - static_cast<double>(n) * result.z) / argument;
        }
        return result;
    };

    for (size_t point = 0; point < x.size(); ++point) {
        const double radius = std::sqrt(x[point] * x[point] + y[point] * y[point] + z[point] * z[point]);
        if (radius < 1e-15) {
            if (!scattered_only && field_type == "Ex") values[point] = e0x;
            else if (!scattered_only && field_type == "Ey") values[point] = e0y;
            continue;
        }
        if (scattered_only && radius < outer_radius) continue;

        const double theta = std::acos(clamp(z[point] / radius));
        const double phi = std::atan2(y[point], x[point]);
        const double st = std::sin(theta), ct = std::cos(theta);
        const double sp = std::sin(phi), cp = std::cos(phi);
        const bool in_core = radius < core_radius;
        const bool in_shell = radius >= core_radius && radius < outer_radius;
        const C argument = in_core
            ? source->wavenumber_vacuum * core_index * radius
            : source->wavenumber_vacuum * (in_shell ? shell_index : scatterer.medium->get_refractive_index()) * radius;
        auto [pi, tau] = scatterer.get_pi_tau(std::cos(theta), scatterer.max_order);

        C erx{}, etx{}, epx{}, ery{}, ety{}, epy{};
        for (size_t index = 0; index < scatterer.max_order; ++index) {
            const size_t n = index + 1;
            const double nd = static_cast<double>(n);
            const C en = C{std::cos(Constants::PI * n / 2.0), std::sin(Constants::PI * n / 2.0)}
                * (2.0 * nd + 1.0) / (nd * (nd + 1.0));
            const C nn1 = nd * (nd + 1.0);
            const C p1 = pi[index] * st;
            C mx, nx, my, ny, md, ndv;
            if (scattered_only) {
                const auto f = radial(n, C{k_medium * radius, 0.0}, 2);
                mx = -scatterer.bn[index] * f.z; nx = -scatterer.an[index] * f.z;
                md = -scatterer.bn[index] * f.derivative; ndv = -scatterer.an[index] * f.derivative;
                my = mx; ny = nx;
            } else if (in_core) {
                const auto f = radial(n, argument, 0);
                mx = scatterer.cn[index] * f.z; nx = scatterer.dn[index] * f.z;
                md = scatterer.cn[index] * f.derivative; ndv = scatterer.dn[index] * f.derivative;
                my = mx; ny = nx;
            } else if (in_shell) {
                const auto f = radial(n, argument, 0);
                const auto fy = radial(n, argument, 1);
                mx = scatterer.shell_m_regular[index] * f.z + scatterer.shell_m_irregular[index] * fy.z;
                nx = scatterer.shell_n_regular[index] * f.z + scatterer.shell_n_irregular[index] * fy.z;
                md = scatterer.shell_m_regular[index] * f.derivative + scatterer.shell_m_irregular[index] * fy.derivative;
                ndv = scatterer.shell_n_regular[index] * f.derivative + scatterer.shell_n_irregular[index] * fy.derivative;
                my = mx; ny = nx;
            } else {
                const auto fj = radial(n, C{k_medium * radius, 0.0}, 0);
                const auto fh = radial(n, C{k_medium * radius, 0.0}, 2);
                mx = fj.z - scatterer.bn[index] * fh.z; nx = fj.z - scatterer.an[index] * fh.z;
                md = fj.derivative - scatterer.bn[index] * fh.derivative;
                ndv = fj.derivative - scatterer.an[index] * fh.derivative;
                my = mx; ny = nx;
            }
            const C mor_t = cp * pi[index] * mx, mor_p = -sp * tau[index] * mx;
            const C mer_t = sp * pi[index] * my, mer_p = cp * tau[index] * my;
            const C ner_r = cp * nn1 * p1 * (nx / argument);
            const C ner_t = cp * tau[index] * ndv, ner_p = -sp * pi[index] * ndv;
            const C nor_r = sp * nn1 * p1 * (ny / argument);
            const C nor_t = sp * tau[index] * ndv, nor_p = cp * pi[index] * ndv;
            erx += en * (-i_unit * ner_r); etx += en * (mor_t - i_unit * ner_t); epx += en * (mor_p - i_unit * ner_p);
            ery += en * (-i_unit * nor_r); ety += en * (mer_t - i_unit * nor_t); epy += en * (mer_p - i_unit * nor_p);
        }
        C er = e0x * erx + e0y * ery, et = e0x * etx + e0y * ety, ep = e0x * epx + e0y * epy;
        const C ex = er * st * cp + et * ct * cp - ep * sp;
        const C ey = er * st * sp + et * ct * sp + ep * cp;
        const C ez = er * ct - et * st;
        if (field_type == "Ex") values[point] = ex;
        else if (field_type == "Ey") values[point] = ey;
        else if (field_type == "Ez") values[point] = ez;
        else values[point] = std::sqrt(std::norm(ex) + std::norm(ey) + std::norm(ez));
    }
    return values;
}
}

std::vector<complex128> CoreShell::get_total_nearfields(
    const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
    const std::string& field_type, const std::shared_ptr<BaseSource>& source) {
    this->compute_cn_dn(this->max_order);
    return evaluate_core_shell_field(*this, x, y, z, field_type, source, false);
}

std::vector<complex128> CoreShell::get_scattered_nearfields(
    const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z,
    const std::string& field_type, const std::shared_ptr<BaseSource>& source) {
    this->compute_an_bn(this->max_order);
    return evaluate_core_shell_field(*this, x, y, z, field_type, source, true);
}
