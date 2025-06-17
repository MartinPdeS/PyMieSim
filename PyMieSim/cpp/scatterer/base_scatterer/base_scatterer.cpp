#include "scatterer/base_scatterer/base_scatterer.h"



size_t BaseScatterer::get_wiscombe_criterion(const double size_parameter) const {
    return static_cast<size_t>(2 + size_parameter + 4 * std::cbrt(size_parameter)) + 16;
}

std::vector<double> BaseScatterer::get_prefactor() const {
    std::vector<double> output;

    output.reserve(max_order);

    for (size_t m = 0; m < max_order ; ++m)
        output[m] = (double) ( 2 * (m + 1) + 1 ) / ( (m + 1) * ( (m + 1) + 1 ) );

    return output;
}


std::tuple<std::vector<complex128>, std::vector<complex128>>
BaseScatterer::compute_structured_fields(const std::vector<complex128>& S1, const std::vector<complex128>& S2, const std::vector<double>& theta, const double& radius) const {
    std::vector<complex128> phi_field, theta_field;

    size_t full_size = theta.size() * S1.size();

    phi_field.reserve(full_size);
    theta_field.reserve(full_size);

    complex128 propagator = this->get_propagator(radius);

    for (unsigned int p=0; p < S1.size(); p++ )
        for (unsigned int t=0; t < theta.size(); t++ )
        {
            complex128 phi_point_field = propagator * S1[p] * (source.jones_vector[0] * cos(theta[t]) + source.jones_vector[1] * sin(theta[t]));
            complex128 thetea_point_field = propagator * S2[p] * (source.jones_vector[0] * sin(theta[t]) - source.jones_vector[1] * cos(theta[t]));

            phi_field.push_back(phi_point_field);
            theta_field.push_back(thetea_point_field);
        }

    return std::make_tuple(phi_field, theta_field);
}

std::tuple<std::vector<complex128>, std::vector<complex128>>
BaseScatterer::compute_structured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const {
    auto [S1, S2] = this->compute_s1s2(phi);

    return this->compute_structured_fields(S1, S2, theta, radius);
}

std::tuple<std::vector<complex128>, std::vector<complex128>>
BaseScatterer::compute_unstructured_fields(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const
{
    auto [S1, S2] = this->compute_s1s2(phi);

    std::vector<complex128> phi_field, theta_field;

    size_t full_size = theta.size();

    phi_field.reserve(full_size);
    theta_field.reserve(full_size);

    complex128 propagator = this->get_propagator(radius);

    for (unsigned int idx=0; idx < full_size; idx++)
    {
        complex128 phi_field_point = propagator * S1[idx] * (source.jones_vector[0] * cos(theta[idx]) + source.jones_vector[1] * sin(theta[idx])),
        theta_field_point = propagator * S2[idx] * (source.jones_vector[0] * sin(theta[idx]) - source.jones_vector[1] * cos(theta[idx]));

        phi_field.push_back(phi_field_point);
        theta_field.push_back(theta_field_point);
    }

    return std::make_tuple(phi_field, theta_field);
}

std::tuple<std::vector<complex128>, std::vector<complex128>>
BaseScatterer::compute_unstructured_fields(const FibonacciMesh& fibonacci_mesh, const double radius) const
{
    return this->compute_unstructured_fields(
        fibonacci_mesh.spherical.phi,
        fibonacci_mesh.spherical.theta,
        radius
    );
}

std::tuple<std::vector<complex128>, std::vector<complex128>, std::vector<double>, std::vector<double>>
BaseScatterer::compute_full_structured_fields(const size_t& sampling, const double& radius) const
{
    FullSteradian full_mesh = FullSteradian(sampling, radius);

    auto [S1, S2] = this->compute_s1s2(full_mesh.spherical.phi);

    auto [phi_field, theta_field] = this->compute_structured_fields(S1, S2, full_mesh.spherical.theta, radius);

    return std::make_tuple(
        std::move(phi_field),
        std::move(theta_field),
        std::move(full_mesh.spherical.phi),
        std::move(full_mesh.spherical.theta)
    );
}

std::vector<complex128>
BaseScatterer::compute_dn(double nmx, complex128 z)  const { //Page 127 of BH
    std::vector<complex128> Dn(nmx, 0.0);

    for (double n = nmx - 1.; n > 1.; n--)
        Dn[n-1] = n/z - ( 1. / (Dn[n] + n/z) );

    return Dn;
}

double
BaseScatterer::get_g_with_fields(size_t sampling) const {
    auto [SPF, fibonacci_mesh] = this->compute_full_structured_spf(sampling);

    double
    norm = abs(fibonacci_mesh.get_integral(SPF)),
    expected_cos = abs(fibonacci_mesh.get_cos_integral(SPF));

    return expected_cos/norm;
}

//- PYTHON INTERFACE --------------------------------------------------------------------------

std::tuple<py::array_t<complex128>, py::array_t<complex128>>
BaseScatterer::get_unstructured_fields_py(const std::vector<double>& phi, const std::vector<double>& theta, const double radius) const
{
    auto [theta_field, phi_field] = this->compute_unstructured_fields(phi, theta, radius);

    return std::make_tuple(
        _vector_to_numpy(phi_field, {phi_field.size()}),
        _vector_to_numpy(theta_field, {theta_field.size()})
    );
}

std::tuple<py::array_t<complex128>, py::array_t<complex128>>
BaseScatterer::get_s1s2_py(const std::vector<double> &phi) const
{
    auto [S1, S2] = this->compute_s1s2(phi);

    return std::make_tuple(
        _vector_to_numpy(S1, {S1.size()}),
        _vector_to_numpy(S2, {S2.size()})
    );
}

std::tuple<py::array_t<complex128>, py::array_t<complex128>, py::array_t<double>, py::array_t<double>>
BaseScatterer::get_full_structured_fields_py(size_t &sampling, double& distance) const {
    auto [phi_field, theta_field, theta, phi] = this->compute_full_structured_fields(sampling, distance);

    py::array_t<complex128>
        phi_field_py = _vector_to_numpy(phi_field, {sampling, sampling}),
        theta_field_py = _vector_to_numpy(theta_field, {sampling, sampling});

    py::array_t<double>
        theta_py = _vector_to_numpy(theta, {theta.size()}),
        phi_py = _vector_to_numpy(phi, {phi.size()});

    phi_field_py = phi_field_py.attr("transpose")();
    theta_field_py = theta_field_py.attr("transpose")();

    return std::make_tuple(phi_field_py, theta_field_py, phi_py, theta_py);
}

std::tuple<std::vector<double>, FullSteradian>
BaseScatterer::compute_full_structured_spf(const size_t sampling, const double radius) const
{
    FullSteradian full_mesh = FullSteradian(sampling, radius);

    auto [phi_field, theta_field] = this->compute_structured_fields(
        full_mesh.spherical.phi,
        full_mesh.spherical.theta,
        radius
    );

    std::vector<double> spf;
    spf.reserve(phi_field.size());

    for (size_t iter = 0; iter < phi_field.size(); iter++)
    {
        double value = pow(abs(phi_field[iter]), 2) + pow(abs(theta_field[iter]), 2);
        spf.push_back( value );
    }

    return std::make_tuple(std::move(spf), std::move(full_mesh));
}

complex128
BaseScatterer::get_coefficient_py(const std::string &type, const size_t order) {
    if (order > max_order)
        throw std::invalid_argument("Coefficient number is higher than computed max value: " + std::to_string(max_order));

    if (type == "a")
        return this->an[order];
    else if (type == "b")
        return this->bn[order];
    else if (type == "c")
        return this->cn[order];
    else if (type == "d")
        return this->dn[order];

    throw std::invalid_argument("Invalid type provided.");

    return 0.;
}

complex128
BaseScatterer::get_propagator(const double &radius) const {
    return source.amplitude / (source.wavenumber * radius) * exp(-complex128(0, 1) * source.wavenumber * radius);
}

std::tuple<std::vector<complex128>, std::vector<complex128>>
BaseScatterer::get_pi_tau(const double& mu, const size_t& max_order) const {
    std::vector<complex128> pin, taun;
    pin.reserve(max_order);
    taun.reserve(max_order);

    pin.push_back( 1. );
    pin.push_back( 3. * mu );

    taun.push_back( mu );
    taun.push_back( 3.0 * cos(2. * acos(mu) ) );

    for (size_t order = 2; order < max_order; order++) {
        pin.push_back( ( (2. * (double)order + 1.) * mu * pin[order - 1] - ((double)order + 1.) * pin[order - 2] ) / (double)order );

        taun.push_back( ((double)order + 1.) * mu * pin[order] - ((double)order + 2.) * pin[order - 1] );
    }

    return std::make_tuple(pin, taun);
}

void
BaseScatterer::get_pi_tau(double mu, size_t max_order, complex128 *pin, complex128 *taun) const {
    pin[0] = 1.;
    pin[1] = 3. * mu;

    taun[0] = mu;
    taun[1] = 3.0 * cos(2. * acos(mu) );

    for (size_t order = 2; order < max_order; order++) {
        pin[order] = ( (2. * (double)order + 1.) * mu * pin[order - 1] - ((double)order + 1.) * pin[order - 2] ) / (double)order;

        taun[order] = ((double)order + 1.) * mu * pin[order] - ((double)order + 2.) * pin[order - 1];
    }
}
