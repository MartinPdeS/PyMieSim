#include "base_scatterer.h"



size_t BaseScatterer::get_wiscombe_criterion(const double size_parameter) const {
    return static_cast<size_t>(2 + size_parameter + 4 * std::cbrt(size_parameter)) + 16;
}

std::vector<double> BaseScatterer::get_prefactor() const {
    std::vector<double> output;
    output.reserve(max_order);

    for (size_t m = 0; m < max_order; ++m) {
        const double n = static_cast<double>(m + 1);
        output.push_back((2.0 * n + 1.0) / (n * (n + 1.0)));
    }

    return output;
}


std::tuple<
    std::vector<complex128>,
    std::vector<complex128>,
    FullSteradian
>
BaseScatterer::get_structured_farfields(
    const size_t& sampling,
    const double& radius,
    const std::shared_ptr<BaseSource>& source
) const {

    FullSteradian mesh(sampling);

    std::vector<complex128> phi_field, theta_field;
    phi_field.reserve(mesh.total_size);
    theta_field.reserve(mesh.total_size);

    std::pair<std::vector<complex128>, std::vector<complex128>> s1s2 = this->compute_s1s2(mesh.spherical.phi);

    complex128 propagator = this->get_propagator(radius, source);

    complex128 E0x = source->polarization.jones_vector[0];
    complex128 E0y = source->polarization.jones_vector[1];

    const std::vector<complex128>& S1 = s1s2.first;
    const std::vector<complex128>& S2 = s1s2.second;

    for (unsigned int p=0; p < sampling; p++ )
        for (unsigned int t=0; t < sampling; t++ )
        {
            complex128 phi_point_field = propagator * S1[p] * (E0x * cos(mesh.spherical.theta[t]) + E0y * sin(mesh.spherical.theta[t]));
            complex128 theta_point_field = propagator * S2[p] * (E0x * sin(mesh.spherical.theta[t]) - E0y * cos(mesh.spherical.theta[t]));

            phi_field.push_back(phi_point_field);
            theta_field.push_back(theta_point_field);
        }

    return std::make_tuple(phi_field, theta_field, mesh);
}

std::pair<std::vector<complex128>, std::vector<complex128>>
BaseScatterer::get_unstructured_farfields(
    const std::vector<double>& phi,
    const std::vector<double>& theta,
    const double radius,
    const std::shared_ptr<BaseSource>& source
) const
{
    auto [S1, S2] = this->compute_s1s2(phi);

    std::vector<complex128> phi_field, theta_field;

    size_t full_size = theta.size();

    phi_field.reserve(full_size);
    theta_field.reserve(full_size);

    complex128 propagator = this->get_propagator(radius, source);

    complex128 E0x = source->polarization.jones_vector[0];
    complex128 E0y = source->polarization.jones_vector[1];

    for (unsigned int idx=0; idx < full_size; idx++)
    {
        complex128 phi_field_point = propagator * S1[idx] * (E0x * cos(theta[idx]) + E0y * sin(theta[idx]));
        complex128 theta_field_point = propagator * S2[idx] * (E0x * sin(theta[idx]) - E0y * cos(theta[idx]));

        phi_field.push_back(phi_field_point);
        theta_field.push_back(theta_field_point);
    }

    return std::make_pair(phi_field, theta_field);
}

std::pair<std::vector<complex128>, std::vector<complex128>>
BaseScatterer::get_unstructured_farfields(
    const FibonacciMesh& fibonacci_mesh,
    const double radius,
    const std::shared_ptr<BaseSource>& source
) const
{
    return this->get_unstructured_farfields(
        fibonacci_mesh.spherical.phi,
        fibonacci_mesh.spherical.theta,
        radius,
        source
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
BaseScatterer::get_g_with_farfields(std::shared_ptr<BaseSource> source, size_t sampling) const {
    auto [SPF, mesh] = this->get_structured_spf(source, sampling);

    double
    norm = abs(mesh.get_integral(SPF)),
    expected_cos = abs(mesh.get_cos_integral(SPF));

    return expected_cos / norm;
}

complex128
BaseScatterer::get_coefficient(const std::string &type, const size_t order) {
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
BaseScatterer::get_propagator(const double &radius, std::shared_ptr<BaseSource> source) const {
    return (
        source->amplitude /
        (source->wavenumber_vacuum * this->medium->get_refractive_index(source->wavelength) * radius) *
        exp(-complex128(0, 1) * source->wavenumber_vacuum * radius)
    );
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

std::pair<std::vector<complex128>, std::vector<complex128>>
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

    return std::make_pair(pin, taun);
}

std::pair<std::vector<double>, FullSteradian>
BaseScatterer::get_structured_spf(std::shared_ptr<BaseSource> source, const size_t sampling, const double radius) const
{
    auto [phi_field, theta_field, full_mesh] = this->get_structured_farfields(
        sampling,
        radius,
        source
    );

    std::vector<double> spf;
    spf.reserve(phi_field.size());

    for (size_t iter = 0; iter < phi_field.size(); iter++)
    {
        double value = pow(abs(phi_field[iter]), 2) + pow(abs(theta_field[iter]), 2);
        spf.push_back( value );
    }

    return std::make_pair(std::move(spf), std::move(full_mesh));
}

std::vector<double> BaseScatterer::get_unstructured_spf(
    std::shared_ptr<BaseSource> source,
    const std::vector<double>& phi,
    const std::vector<double>& theta,
    const double radius
) const {
    auto [phi_field, theta_field] = this->get_unstructured_farfields(
        phi,
        theta,
        radius,
        source
    );

    std::vector<double> spf;
    spf.reserve(phi_field.size());

    for (size_t iter = 0; iter < phi_field.size(); iter++)
    {
        double value = pow(abs(phi_field[iter]), 2) + pow(abs(theta_field[iter]), 2);
        spf.push_back( value );
    }

    return spf;
}

std::tuple<
    std::vector<complex128>,
    std::vector<complex128>,
    std::vector<complex128>,
    std::vector<complex128>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
>
BaseScatterer::get_structured_total_nearfields(
    const std::vector<double>& x_range,
    const std::vector<double>& y_range,
    const std::vector<double>& z_range,
    const std::string& field_type,
    std::shared_ptr<BaseSource> source
) {
    const size_t number_of_x_points = x_range.size();
    const size_t number_of_y_points = y_range.size();
    const size_t number_of_z_points = z_range.size();

    if (number_of_x_points == 0 || number_of_y_points == 0 || number_of_z_points == 0) {
        throw std::runtime_error("x_range, y_range, and z_range must all be non empty.");
    }

    if (number_of_x_points > std::numeric_limits<size_t>::max() / number_of_y_points) {
        throw std::runtime_error("Structured grid size overflow while computing number_of_x_points * number_of_y_points.");
    }

    const size_t number_of_xy_points = number_of_x_points * number_of_y_points;

    if (number_of_xy_points > std::numeric_limits<size_t>::max() / number_of_z_points) {
        throw std::runtime_error("Structured grid size overflow while computing total number of grid points.");
    }

    const size_t total_number_of_points = number_of_xy_points * number_of_z_points;

    std::vector<double> x_coords(total_number_of_points);
    std::vector<double> y_coords(total_number_of_points);
    std::vector<double> z_coords(total_number_of_points);

    #pragma omp parallel for collapse(3)
    for (size_t x_index = 0; x_index < number_of_x_points; ++x_index) {
        for (size_t y_index = 0; y_index < number_of_y_points; ++y_index) {
            for (size_t z_index = 0; z_index < number_of_z_points; ++z_index) {
                const size_t linear_index =
                    x_index * (number_of_y_points * number_of_z_points) +
                    y_index * number_of_z_points +
                    z_index;

                x_coords[linear_index] = x_range[x_index];
                y_coords[linear_index] = y_range[y_index];
                z_coords[linear_index] = z_range[z_index];
            }
        }
    }

    std::vector<complex128> field_values =
        this->get_total_nearfields(x_coords, y_coords, z_coords, field_type, source);

    std::vector<complex128> field_x_components =
        this->get_total_nearfields(x_coords, y_coords, z_coords, "Ex", source);

    std::vector<complex128> field_y_components =
        this->get_total_nearfields(x_coords, y_coords, z_coords, "Ey", source);

    std::vector<complex128> field_z_components =
        this->get_total_nearfields(x_coords, y_coords, z_coords, "Ez", source);

    return std::make_tuple(
        std::move(field_values),
        std::move(field_x_components),
        std::move(field_y_components),
        std::move(field_z_components),
        std::move(x_coords),
        std::move(y_coords),
        std::move(z_coords)
    );
}

std::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>
>
BaseScatterer::get_unstructured_stokes(
    const std::vector<double>& phi,
    const std::vector<double>& theta,
    const double distance,
    std::shared_ptr<BaseSource> source
) const
{
    auto [E_phi, E_theta] = this->get_unstructured_farfields(phi, theta, distance, source);

    std::vector<double> I(E_phi.size());
    std::vector<double> Q(E_phi.size());
    std::vector<double> U(E_phi.size());
    std::vector<double> V(E_phi.size());

    double I_max = 0.0;

    for (size_t i = 0; i < E_phi.size(); ++i) {
        const double intensity = std::norm(E_phi[i]) + std::norm(E_theta[i]);
        I[i] = intensity;
        I_max = std::max(I_max, intensity);

        Q[i] = (std::norm(E_phi[i]) - std::norm(E_theta[i])) / intensity;
        U[i] = (2.0 * std::real(E_phi[i] * std::conj(E_theta[i]))) / intensity;
        V[i] = (-2.0 * std::imag(E_phi[i] * std::conj(E_theta[i]))) / intensity;
    }

    if (I_max > 0.0) {
        for (double& v : I) v /= I_max;
    }

    return std::make_tuple(
        std::move(I),
        std::move(Q),
        std::move(U),
        std::move(V)
    );
}

std::tuple<
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    std::vector<double>,
    FullSteradian
>
BaseScatterer::get_structured_stokes(
    const size_t sampling,
    const double distance,
    std::shared_ptr<BaseSource> source
) const {
    auto [E_phi, E_theta, full_mesh] = this->get_structured_farfields(sampling, distance, source);

    size_t total_size = full_mesh.sampling * full_mesh.sampling;
    std::vector<double> I(total_size);
    std::vector<double> Q(total_size);
    std::vector<double> U(total_size);
    std::vector<double> V(total_size);

    double I_max = 0.0;

    for (size_t i = 0; i < total_size; ++i) {
        const double intensity = std::norm(E_phi[i]) + std::norm(E_theta[i]);
        I[i] = intensity;
        I_max = std::max(I_max, intensity);

        Q[i] = (std::norm(E_phi[i]) - std::norm(E_theta[i])) / intensity;
        U[i] = (2.0 * std::real(E_phi[i] * std::conj(E_theta[i]))) / intensity;
        V[i] = (-2.0 * std::imag(E_phi[i] * std::conj(E_theta[i]))) / intensity;
    }

    if (I_max > 0.0) {
        for (double& v : I)
            v /= I_max;
    }

    return std::make_tuple(
        std::move(I),
        std::move(Q),
        std::move(U),
        std::move(V),
        std::move(full_mesh)
    );
}


std::vector<complex128>
BaseScatterer::compute_incident_nearfields(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const std::vector<double>& z,
    const std::string& field_type,
    const std::shared_ptr<BaseSource>& source
) const
{
    if (field_type != "Ex" && field_type != "Ey" && field_type != "Ez" && field_type != "|E|") {
        throw std::invalid_argument("Invalid field_type. Must be one of: Ex, Ey, Ez, |E|");
    }

    const std::size_t number_of_points = x.size();
    if (number_of_points != y.size() || number_of_points != z.size()) {
        throw std::invalid_argument("x, y, z vectors must have the same length");
    }

    std::vector<complex128> field_values(number_of_points);

    const double k0 = source->wavenumber_vacuum;
    const double k_medium = k0 * this->medium->get_refractive_index(source->wavelength);

    const complex128 E0x = source->polarization.jones_vector[0] * source->amplitude;
    const complex128 E0y = source->polarization.jones_vector[1] * source->amplitude;

    const complex128 i_unit(0.0, 1.0);

    for (std::size_t point_index = 0; point_index < number_of_points; ++point_index) {

        // const double x_position = x[point_index];
        // const double y_position = y[point_index];
        const double z_position = z[point_index];

        // Plane wave consistent with your BH expansion assumptions: propagation along +z in the medium
        const complex128 phase = std::exp(i_unit * complex128(k_medium * z_position, 0.0));

        const complex128 Ex_inc = E0x * phase;
        const complex128 Ey_inc = E0y * phase;
        const complex128 Ez_inc(0.0, 0.0);

        if (field_type == "Ex") {
            field_values[point_index] = Ex_inc;
        } else if (field_type == "Ey") {
            field_values[point_index] = Ey_inc;
        } else if (field_type == "Ez") {
            field_values[point_index] = Ez_inc;
        } else {
            const double mag = std::sqrt(std::norm(Ex_inc) + std::norm(Ey_inc) + std::norm(Ez_inc));
            field_values[point_index] = complex128(mag, 0.0);
        }
    }

    return field_values;
}
