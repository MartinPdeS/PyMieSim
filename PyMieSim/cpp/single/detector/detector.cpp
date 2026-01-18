#include "detector.h"

// ------------------------- Initialization Function -------------------------
void Photodiode::initialize()
{
    this->parse_mode(this->mode_number);
    this->mode_field = ModeField(this->mode_id);

    this->snell_interface.set_media(
        this->scatterer_medium_refractive_index,
        this->detector_medium_refractive_index
    );

    double theta_s_max = 0.0;
    double theta_s_min = 0.0;

    this->snell_interface.compute_snell_mapped_angle_bounds(
        this->numerical_aperture,
        this->cache_numerical_aperture,
        theta_s_max,
        theta_s_min
    );

    this->max_angle = theta_s_max;
    this->min_angle = theta_s_min;

    this->fibonacci_mesh = FibonacciMesh(
        this->sampling,
        this->max_angle,
        this->min_angle,
        this->phi_offset,
        this->gamma_offset,
        this->rotation
    );

    this->scalar_field = this->mode_field.get_unstructured(
        this->fibonacci_mesh.base_cartesian.x,
        this->fibonacci_mesh.base_cartesian.y
    );

    this->snell_interface.compute_fresnel_amplitude_transmission_on_mesh(
        this->fibonacci_mesh,
        this->interface_t_s,
        this->interface_t_p
    );
}




// ------------------------- Utils Functions -------------------------
void Photodiode::parse_mode(const std::string& mode_number) {
    if (mode_number.size() < 4) {
        throw std::invalid_argument("mode string too short, need at least 4 chars");
    }

    this->mode_id.mode_family = mode_number.substr(0, 2);
    char c0 = mode_number[2], c1 = mode_number[3];
    if (!std::isdigit(c0) || !std::isdigit(c1)) {
        throw std::invalid_argument("expected digits in positions 2 and 3");
    }
    this->mode_id.number_0 = c0 - '0';
    this->mode_id.number_1 = c1 - '0';
}

template <class T> inline
double Photodiode::get_norm1_squared(const std::vector<T> &array) const
{
    T sum  = 0.0;

    for (auto v : array)
        sum += v;

    return pow(abs(sum), 2);
}

template <class T> inline
double Photodiode::get_norm2_squared(const std::vector<T> &array) const
{
    T sum  = 0.0;

    for (auto v : array)
    sum += pow( abs(v), 2 );

    return abs(sum);
}


template <class T> inline
void Photodiode::square_array(std::vector<T> &array)
{
    for (T &v : array)
        v = pow( abs(v), 2);
}


std::tuple<std::vector<complex128>, std::vector<complex128>>
Photodiode::get_projected_farfields(const std::vector<complex128> &theta_field, const std::vector<complex128> &phi_field) const
{

    std::vector<complex128>
        horizontal_projection(theta_field.size()),
        vertical_projection(theta_field.size());

    for (size_t i=0; i<theta_field.size(); ++i)
    {
        vertical_projection[i] =
            theta_field[i] * this->fibonacci_mesh.vertical_perpendicular_projection[i] +
            phi_field[i] * this->fibonacci_mesh.vertical_parallel_projection[i] ;  // new_version



        horizontal_projection[i] =
            theta_field[i] * this->fibonacci_mesh.horizontal_perpendicular_projection[i] +
            phi_field[i] * this->fibonacci_mesh.horizontal_parallel_projection[i] ; // new_version
    }

    return std::make_tuple(horizontal_projection, vertical_projection);

}


void Photodiode::apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const //Theta = Para
{
    for (size_t i=0; i<field0.size(); i++) {
        field0[i] *= this->scalar_field[i];
        field1[i] *= this->scalar_field[i];
    }
}

template <typename T> inline
void Photodiode::apply_polarization_filter(T &coupling_theta, T &coupling_phi, double polarization_filter) const
{

    if (std::isnan(this->polarization_filter))
        return;

    double
        theta_polarization_filtering  = pow( sin(polarization_filter), 2 ),
        phi_polarization_filtering    = pow( cos(polarization_filter), 2 );

    coupling_theta *= theta_polarization_filtering;
    coupling_phi *= phi_polarization_filtering;
}


std::vector<complex128> Photodiode::get_structured_scalarfield(const size_t sampling) const {
    // Define the spatial extent (adjust as needed)
    const double xmin = -1.0;
    const double xmax = 1.0;
    const double ymin = -1.0;
    const double ymax = 1.0;

    std::vector<double> x_coords(sampling);
    std::vector<double> y_coords(sampling);
    const double dx = (xmax - xmin) / (sampling - 1);
    const double dy = (ymax - ymin) / (sampling - 1);

    for (size_t i = 0; i < sampling; ++i) {
        x_coords[i] = xmin + i * dx;
        y_coords[i] = ymin + i * dy;
    }

    // Generate structured meshgrid (flattened)
    std::vector<double> X_flat, Y_flat;
    X_flat.reserve(sampling * sampling);
    Y_flat.reserve(sampling * sampling);

    for (size_t j = 0; j < sampling; ++j) {
        for (size_t i = 0; i < sampling; ++i) {
                X_flat.push_back(x_coords[i]);  // column (x)
                Y_flat.push_back(y_coords[j]);  // row (y)
        }
    }

    // Compute unstructured field
    std::vector<complex128> output = this->mode_field.get_unstructured(X_flat, Y_flat);

    return output;
}

std::vector<double> Photodiode::get_poynting_field(const BaseScatterer& scatterer, double distance) const {
    auto [Ephi, Etheta] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, distance);

    double electric_field_norm = 0.0;
    double magnetic_field_norm = 0.0;

    std::vector<double> poynting(Ephi.size());

    for (size_t i = 0; i < Ephi.size(); ++i) {
        electric_field_norm = std::sqrt(std::pow(std::abs(Ephi[i]), 2) + std::pow(std::abs(Etheta[i]), 2));
        magnetic_field_norm = (electric_field_norm / LIGHT_SPEED);  // units of Tesla
        poynting[i] = EPSILON0 * std::pow(LIGHT_SPEED, 2) * electric_field_norm * magnetic_field_norm;  // units of W/m^2
    }

    return poynting;
}


double Photodiode::get_energy_flow(const BaseScatterer& scatterer) const {
    const double distance = 1;  // Energy flow do not depend on distance
    std::vector<double> poynting_vector = this->get_poynting_field(scatterer, distance);

    double dA = distance * distance * this->fibonacci_mesh.dOmega;  // units of m^2

    return 0.5 * trapz(poynting_vector, dA); // units of W
}

// ------------------------- Coupling Function -------------------------
double Photodiode::get_coupling(const BaseScatterer& scatterer) const {
    if (this->is_coherent)
        return this->mean_coupling ? get_coupling_mean_coherent(scatterer) : get_coupling_point_coherent(scatterer);
    else
        return this->mean_coupling ? get_coupling_mean_no_coherent(scatterer) : get_coupling_point_no_coherent(scatterer);
}


double Photodiode::get_coupling_mean_coherent(const BaseScatterer& scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    // New line
    this->apply_interface_transmission_to_fields(theta_field, phi_field);

    auto [horizontal_projection, vertical_projection] = this->get_projected_farfields(theta_field, phi_field);

    this->apply_scalar_field(horizontal_projection, vertical_projection);

    double coupling_theta = this->get_norm2_squared(horizontal_projection);
    double coupling_phi   = this->get_norm2_squared(vertical_projection);

    this->apply_polarization_filter(coupling_theta, coupling_phi, this->polarization_filter);

    return 0.5 * EPSILON0 * LIGHT_SPEED
         * (coupling_theta + coupling_phi)
         * this->fibonacci_mesh.dOmega / this->fibonacci_mesh.Omega;
}


double Photodiode::get_coupling_mean_no_coherent(const BaseScatterer &scatterer) const
{
    return get_coupling_point_no_coherent(scatterer);
}

double Photodiode::get_coupling_point_coherent(const BaseScatterer& scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    this->apply_interface_transmission_to_fields(theta_field, phi_field);

    auto [horizontal_projection, vertical_projection] = this->get_projected_farfields(theta_field, phi_field);

    this->apply_scalar_field(horizontal_projection, vertical_projection);

    double coupling_theta = this->get_norm1_squared(horizontal_projection);
    double coupling_phi   = this->get_norm1_squared(vertical_projection);

    this->apply_polarization_filter(coupling_theta, coupling_phi, this->polarization_filter);

    return 0.5 * EPSILON0 * LIGHT_SPEED
         * (coupling_theta + coupling_phi)
         * this->fibonacci_mesh.dOmega;
}

double Photodiode::get_coupling_point_no_coherent(const BaseScatterer& scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    this->apply_interface_transmission_to_fields(theta_field, phi_field);

    double coupling_theta = this->get_norm2_squared(theta_field);
    double coupling_phi   = this->get_norm2_squared(phi_field);

    this->apply_polarization_filter(coupling_theta, coupling_phi, this->polarization_filter);

    return 0.5 * EPSILON0 * LIGHT_SPEED
         * (coupling_theta + coupling_phi)
         * this->fibonacci_mesh.dOmega;
}


// detector.cpp (additions)

void Photodiode::apply_interface_transmission_to_fields(
    std::vector<complex128>& theta_field,
    std::vector<complex128>& phi_field
) const
{
    if (this->interface_t_s.size() != theta_field.size() ||
        this->interface_t_p.size() != phi_field.size()) {
        throw std::runtime_error("interface transmission not initialized or size mismatch.");
    }

    // Interface normal (detector axis)
    const double nx = 0.0;
    const double ny = 0.0;
    const double nz = 1.0;

    for (size_t i = 0; i < theta_field.size(); ++i) {

        const double kx = this->fibonacci_mesh.cartesian.x[i];
        const double ky = this->fibonacci_mesh.cartesian.y[i];
        const double kz = this->fibonacci_mesh.cartesian.z[i];

        // Guard: only forward hemisphere
        if (kz <= 0.0) {
            theta_field[i] = complex128(0.0, 0.0);
            phi_field[i]   = complex128(0.0, 0.0);
            continue;
        }

        // Convert direction to spherical angles
        const double theta = std::acos(clamp_m1_p1(kz));
        const double phi   = std::atan2(ky, kx);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double sin_phi   = std::sin(phi);
        const double cos_phi   = std::cos(phi);

        // Spherical basis vectors in cartesian coordinates
        // e_theta = (cosθ cosφ, cosθ sinφ, -sinθ)
        // e_phi   = (-sinφ,     cosφ,      0)
        const double e_theta_x = cos_theta * cos_phi;
        const double e_theta_y = cos_theta * sin_phi;
        const double e_theta_z = -sin_theta;

        const double e_phi_x = -sin_phi;
        const double e_phi_y =  cos_phi;
        const double e_phi_z =  0.0;

        // Field in cartesian: E = Eθ eθ + Eφ eφ
        const complex128 Ex = theta_field[i] * e_theta_x + phi_field[i] * e_phi_x;
        const complex128 Ey = theta_field[i] * e_theta_y + phi_field[i] * e_phi_y;
        const complex128 Ez = theta_field[i] * e_theta_z + phi_field[i] * e_phi_z;

        // s direction: s = normalize(n x k)
        double sx = ny * kz - nz * ky; // 0*kz - 1*ky = -ky
        double sy = nz * kx - nx * kz; // 1*kx - 0*kz = kx
        double sz = nx * ky - ny * kx; // 0

        const double s_norm = std::sqrt(sx * sx + sy * sy + sz * sz);

        // If k parallel to n, polarization basis is degenerate. Choose identity.
        if (s_norm == 0.0) {
            continue;
        }

        sx /= s_norm;
        sy /= s_norm;
        sz /= s_norm;

        // p direction: p = s x k
        const double px = sy * kz - sz * ky;
        const double py = sz * kx - sx * kz;
        const double pz = sx * ky - sy * kx;

        // Decompose E into s and p
        const complex128 E_s = Ex * sx + Ey * sy + Ez * sz;
        const complex128 E_p = Ex * px + Ey * py + Ez * pz;

        // Apply Fresnel amplitude transmissions
        const double t_s = this->interface_t_s[i];
        const double t_p = this->interface_t_p[i];

        const complex128 Ex_t = E_s * (t_s * sx) + E_p * (t_p * px);
        const complex128 Ey_t = E_s * (t_s * sy) + E_p * (t_p * py);
        const complex128 Ez_t = E_s * (t_s * sz) + E_p * (t_p * pz);

        // Project back onto spherical basis
        theta_field[i] = Ex_t * e_theta_x + Ey_t * e_theta_y + Ez_t * e_theta_z;
        phi_field[i]   = Ex_t * e_phi_x   + Ey_t * e_phi_y   + Ez_t * e_phi_z;
    }
}

void CoherentMode::initialize()
{
    this->parse_mode(this->mode_number);
    this->mode_field = ModeField(this->mode_id);

    this->snell_interface.set_media(
        this->scatterer_medium_refractive_index,
        this->detector_medium_refractive_index
    );

    double theta_s_max = 0.0;
    double theta_s_min = 0.0;

    this->snell_interface.compute_snell_mapped_angle_bounds(
        this->numerical_aperture,
        this->cache_numerical_aperture,
        theta_s_max,
        theta_s_min
    );

    this->max_angle = theta_s_max;
    this->min_angle = theta_s_min;

    this->fibonacci_mesh = FibonacciMesh(
        this->sampling,
        this->max_angle,
        this->min_angle,
        this->phi_offset,
        this->gamma_offset,
        this->rotation
    );

    this->scalar_field = this->mode_field.get_unstructured(
        this->fibonacci_mesh.base_cartesian.x,
        this->fibonacci_mesh.base_cartesian.y
    );

    this->snell_interface.compute_fresnel_amplitude_transmission_on_mesh(
        this->fibonacci_mesh,
        this->interface_t_s,
        this->interface_t_p
    );
}




// ------------------------- Utils Functions -------------------------
void CoherentMode::parse_mode(const std::string& mode_number) {
    if (mode_number.size() < 4) {
        throw std::invalid_argument("mode string too short, need at least 4 chars");
    }

    this->mode_id.mode_family = mode_number.substr(0, 2);
    char c0 = mode_number[2], c1 = mode_number[3];
    if (!std::isdigit(c0) || !std::isdigit(c1)) {
        throw std::invalid_argument("expected digits in positions 2 and 3");
    }
    this->mode_id.number_0 = c0 - '0';
    this->mode_id.number_1 = c1 - '0';
}

template <class T> inline
double CoherentMode::get_norm1_squared(const std::vector<T> &array) const
{
    T sum  = 0.0;

    for (auto v : array)
        sum += v;

    return pow(abs(sum), 2);
}

template <class T> inline
double CoherentMode::get_norm2_squared(const std::vector<T> &array) const
{
    T sum  = 0.0;

    for (auto v : array)
    sum += pow( abs(v), 2 );

    return abs(sum);
}


template <class T> inline
void CoherentMode::square_array(std::vector<T> &array)
{
    for (T &v : array)
        v = pow( abs(v), 2);
}


std::tuple<std::vector<complex128>, std::vector<complex128>>
CoherentMode::get_projected_farfields(const std::vector<complex128> &theta_field, const std::vector<complex128> &phi_field) const
{

    std::vector<complex128>
        horizontal_projection(theta_field.size()),
        vertical_projection(theta_field.size());

    for (size_t i=0; i<theta_field.size(); ++i)
    {
        vertical_projection[i] =
            theta_field[i] * this->fibonacci_mesh.vertical_perpendicular_projection[i] +
            phi_field[i] * this->fibonacci_mesh.vertical_parallel_projection[i] ;  // new_version



        horizontal_projection[i] =
            theta_field[i] * this->fibonacci_mesh.horizontal_perpendicular_projection[i] +
            phi_field[i] * this->fibonacci_mesh.horizontal_parallel_projection[i] ; // new_version
    }

    return std::make_tuple(horizontal_projection, vertical_projection);

}


void CoherentMode::apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const //Theta = Para
{
    for (size_t i=0; i<field0.size(); i++) {
        field0[i] *= this->scalar_field[i];
        field1[i] *= this->scalar_field[i];
    }
}

template <typename T> inline
void CoherentMode::apply_polarization_filter(T &coupling_theta, T &coupling_phi, double polarization_filter) const
{

    if (std::isnan(this->polarization_filter))
        return;

    double
        theta_polarization_filtering  = pow( sin(polarization_filter), 2 ),
        phi_polarization_filtering    = pow( cos(polarization_filter), 2 );

    coupling_theta *= theta_polarization_filtering;
    coupling_phi *= phi_polarization_filtering;
}


std::vector<complex128> CoherentMode::get_structured_scalarfield(const size_t sampling) const {
    // Define the spatial extent (adjust as needed)
    const double xmin = -1.0;
    const double xmax = 1.0;
    const double ymin = -1.0;
    const double ymax = 1.0;

    std::vector<double> x_coords(sampling);
    std::vector<double> y_coords(sampling);
    const double dx = (xmax - xmin) / (sampling - 1);
    const double dy = (ymax - ymin) / (sampling - 1);

    for (size_t i = 0; i < sampling; ++i) {
        x_coords[i] = xmin + i * dx;
        y_coords[i] = ymin + i * dy;
    }

    // Generate structured meshgrid (flattened)
    std::vector<double> X_flat, Y_flat;
    X_flat.reserve(sampling * sampling);
    Y_flat.reserve(sampling * sampling);

    for (size_t j = 0; j < sampling; ++j) {
        for (size_t i = 0; i < sampling; ++i) {
                X_flat.push_back(x_coords[i]);  // column (x)
                Y_flat.push_back(y_coords[j]);  // row (y)
        }
    }

    // Compute unstructured field
    std::vector<complex128> output = this->mode_field.get_unstructured(X_flat, Y_flat);

    return output;
}

std::vector<double> CoherentMode::get_poynting_field(const BaseScatterer& scatterer, double distance) const {
    auto [Ephi, Etheta] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, distance);

    double electric_field_norm = 0.0;
    double magnetic_field_norm = 0.0;

    std::vector<double> poynting(Ephi.size());

    for (size_t i = 0; i < Ephi.size(); ++i) {
        electric_field_norm = std::sqrt(std::pow(std::abs(Ephi[i]), 2) + std::pow(std::abs(Etheta[i]), 2));
        magnetic_field_norm = (electric_field_norm / LIGHT_SPEED);  // units of Tesla
        poynting[i] = EPSILON0 * std::pow(LIGHT_SPEED, 2) * electric_field_norm * magnetic_field_norm;  // units of W/m^2
    }

    return poynting;
}


double CoherentMode::get_energy_flow(const BaseScatterer& scatterer) const {
    const double distance = 1;  // Energy flow do not depend on distance
    std::vector<double> poynting_vector = this->get_poynting_field(scatterer, distance);

    double dA = distance * distance * this->fibonacci_mesh.dOmega;  // units of m^2

    return 0.5 * trapz(poynting_vector, dA); // units of W
}

// ------------------------- Coupling Function -------------------------
double CoherentMode::get_coupling(const BaseScatterer& scatterer) const {
    if (this->is_coherent)
        return this->mean_coupling ? get_coupling_mean_coherent(scatterer) : get_coupling_point_coherent(scatterer);
    else
        return this->mean_coupling ? get_coupling_mean_no_coherent(scatterer) : get_coupling_point_no_coherent(scatterer);
}


double CoherentMode::get_coupling_mean_coherent(const BaseScatterer& scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    // New line
    this->apply_interface_transmission_to_fields(theta_field, phi_field);

    auto [horizontal_projection, vertical_projection] = this->get_projected_farfields(theta_field, phi_field);

    this->apply_scalar_field(horizontal_projection, vertical_projection);

    double coupling_theta = this->get_norm2_squared(horizontal_projection);
    double coupling_phi   = this->get_norm2_squared(vertical_projection);

    this->apply_polarization_filter(coupling_theta, coupling_phi, this->polarization_filter);

    return 0.5 * EPSILON0 * LIGHT_SPEED
         * (coupling_theta + coupling_phi)
         * this->fibonacci_mesh.dOmega / this->fibonacci_mesh.Omega;
}


double CoherentMode::get_coupling_mean_no_coherent(const BaseScatterer &scatterer) const
{
    return get_coupling_point_no_coherent(scatterer);
}

double CoherentMode::get_coupling_point_coherent(const BaseScatterer& scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    this->apply_interface_transmission_to_fields(theta_field, phi_field);

    auto [horizontal_projection, vertical_projection] = this->get_projected_farfields(theta_field, phi_field);

    this->apply_scalar_field(horizontal_projection, vertical_projection);

    double coupling_theta = this->get_norm1_squared(horizontal_projection);
    double coupling_phi   = this->get_norm1_squared(vertical_projection);

    this->apply_polarization_filter(coupling_theta, coupling_phi, this->polarization_filter);

    return 0.5 * EPSILON0 * LIGHT_SPEED
         * (coupling_theta + coupling_phi)
         * this->fibonacci_mesh.dOmega;
}

double CoherentMode::get_coupling_point_no_coherent(const BaseScatterer& scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    this->apply_interface_transmission_to_fields(theta_field, phi_field);

    double coupling_theta = this->get_norm2_squared(theta_field);
    double coupling_phi   = this->get_norm2_squared(phi_field);

    this->apply_polarization_filter(coupling_theta, coupling_phi, this->polarization_filter);

    return 0.5 * EPSILON0 * LIGHT_SPEED
         * (coupling_theta + coupling_phi)
         * this->fibonacci_mesh.dOmega;
}


// detector.cpp (additions)

void CoherentMode::apply_interface_transmission_to_fields(
    std::vector<complex128>& theta_field,
    std::vector<complex128>& phi_field
) const
{
    if (this->interface_t_s.size() != theta_field.size() ||
        this->interface_t_p.size() != phi_field.size()) {
        throw std::runtime_error("interface transmission not initialized or size mismatch.");
    }

    // Interface normal (detector axis)
    const double nx = 0.0;
    const double ny = 0.0;
    const double nz = 1.0;

    for (size_t i = 0; i < theta_field.size(); ++i) {

        const double kx = this->fibonacci_mesh.cartesian.x[i];
        const double ky = this->fibonacci_mesh.cartesian.y[i];
        const double kz = this->fibonacci_mesh.cartesian.z[i];

        // Guard: only forward hemisphere
        if (kz <= 0.0) {
            theta_field[i] = complex128(0.0, 0.0);
            phi_field[i]   = complex128(0.0, 0.0);
            continue;
        }

        // Convert direction to spherical angles
        const double theta = std::acos(clamp_m1_p1(kz));
        const double phi   = std::atan2(ky, kx);

        const double sin_theta = std::sin(theta);
        const double cos_theta = std::cos(theta);
        const double sin_phi   = std::sin(phi);
        const double cos_phi   = std::cos(phi);

        // Spherical basis vectors in cartesian coordinates
        // e_theta = (cosθ cosφ, cosθ sinφ, -sinθ)
        // e_phi   = (-sinφ,     cosφ,      0)
        const double e_theta_x = cos_theta * cos_phi;
        const double e_theta_y = cos_theta * sin_phi;
        const double e_theta_z = -sin_theta;

        const double e_phi_x = -sin_phi;
        const double e_phi_y =  cos_phi;
        const double e_phi_z =  0.0;

        // Field in cartesian: E = Eθ eθ + Eφ eφ
        const complex128 Ex = theta_field[i] * e_theta_x + phi_field[i] * e_phi_x;
        const complex128 Ey = theta_field[i] * e_theta_y + phi_field[i] * e_phi_y;
        const complex128 Ez = theta_field[i] * e_theta_z + phi_field[i] * e_phi_z;

        // s direction: s = normalize(n x k)
        double sx = ny * kz - nz * ky; // 0*kz - 1*ky = -ky
        double sy = nz * kx - nx * kz; // 1*kx - 0*kz = kx
        double sz = nx * ky - ny * kx; // 0

        const double s_norm = std::sqrt(sx * sx + sy * sy + sz * sz);

        // If k parallel to n, polarization basis is degenerate. Choose identity.
        if (s_norm == 0.0) {
            continue;
        }

        sx /= s_norm;
        sy /= s_norm;
        sz /= s_norm;

        // p direction: p = s x k
        const double px = sy * kz - sz * ky;
        const double py = sz * kx - sx * kz;
        const double pz = sx * ky - sy * kx;

        // Decompose E into s and p
        const complex128 E_s = Ex * sx + Ey * sy + Ez * sz;
        const complex128 E_p = Ex * px + Ey * py + Ez * pz;

        // Apply Fresnel amplitude transmissions
        const double t_s = this->interface_t_s[i];
        const double t_p = this->interface_t_p[i];

        const complex128 Ex_t = E_s * (t_s * sx) + E_p * (t_p * px);
        const complex128 Ey_t = E_s * (t_s * sy) + E_p * (t_p * py);
        const complex128 Ez_t = E_s * (t_s * sz) + E_p * (t_p * pz);

        // Project back onto spherical basis
        theta_field[i] = Ex_t * e_theta_x + Ey_t * e_theta_y + Ez_t * e_theta_z;
        phi_field[i]   = Ex_t * e_phi_x   + Ey_t * e_phi_y   + Ez_t * e_phi_z;
    }
}

