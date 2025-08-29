#include "detector/detector.h"


// ------------------------- Initialization Function -------------------------
void Detector::initialize(const double &medium_refractive_index) {
    this->parse_mode(this->mode_number);
    this->mode_field = ModeField(this->mode_id);

    this->max_angle = this->numercical_aperture_to_angle(this->numerical_aperture / medium_refractive_index);
    this->min_angle = this->numercical_aperture_to_angle(this->cache_numerical_aperture / medium_refractive_index);

    if (this->max_angle < this->min_angle)
        throw std::invalid_argument("Cache NA cannot be larger than detector NA.");


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

}

// ------------------------- Utils Functions -------------------------
void Detector::parse_mode(const std::string& mode_number) {
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


double Detector::numercical_aperture_to_angle(double numerical_aperture) const {
    if (numerical_aperture <= 1.0)
        return asin(numerical_aperture);

    if (numerical_aperture >= 1.0)
        return asin(numerical_aperture - 1.0) + PI / 2.0;

    return 1.0;
}

template <class T> inline
double Detector::get_norm1_squared(const std::vector<T> &array) const
{
    T sum  = 0.0;

    for (auto v : array)
        sum += v;

    return pow(abs(sum), 2);
}

template <class T> inline
double Detector::get_norm2_squared(const std::vector<T> &array) const
{
    T sum  = 0.0;

    for (auto v : array)
    sum += pow( abs(v), 2 );

    return abs(sum);
}


template <class T> inline
void Detector::square_array(std::vector<T> &array)
{
    for (T &v : array)
        v = pow( abs(v), 2);
}


std::tuple<std::vector<complex128>, std::vector<complex128>>
Detector::get_projected_farfields(const std::vector<complex128> &theta_field, const std::vector<complex128> &phi_field) const
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


void Detector::apply_scalar_field(std::vector<complex128> &field0, std::vector<complex128> &field1) const //Theta = Para
{
    for (size_t i=0; i<field0.size(); i++) {
        field0[i] *= this->scalar_field[i];
        field1[i] *= this->scalar_field[i];
    }
}

template <typename T> inline
void Detector::apply_polarization_filter(T &coupling_theta, T &coupling_phi, double polarization_filter) const
{

    if (std::isnan(this->polarization_filter))
        return;

    double
        theta_polarization_filtering  = pow( sin(polarization_filter), 2 ),
        phi_polarization_filtering    = pow( cos(polarization_filter), 2 );

    coupling_theta *= theta_polarization_filtering;
    coupling_phi *= phi_polarization_filtering;
}


pybind11::array_t<complex128> Detector::get_structured_scalarfield(const size_t sampling) const {
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

    // Reshape to 2D numpy array
    return _vector_to_numpy(output, {sampling, sampling});
}


// ------------------------- Coupling Function -------------------------
double Detector::get_coupling(const BaseScatterer& scatterer) const {
    if (this->coherent)
        return this->mean_coupling ? get_coupling_mean_coherent(scatterer) : get_coupling_point_coherent(scatterer);
    else
        return this->mean_coupling ? get_coupling_mean_no_coherent(scatterer) : get_coupling_point_no_coherent(scatterer);
}


double Detector::get_coupling_mean_coherent(const BaseScatterer &scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    auto [horizontal_projection, vertical_projection] = this->get_projected_farfields(theta_field, phi_field);

    this->apply_scalar_field(horizontal_projection, vertical_projection);

    double
        coupling_theta = this->get_norm2_squared(horizontal_projection),
        coupling_phi = this->get_norm2_squared(vertical_projection);

    this->apply_polarization_filter(
        coupling_theta,
        coupling_phi,
        this->polarization_filter
    );

    return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->fibonacci_mesh.dOmega / this->fibonacci_mesh.Omega;
}


double Detector::get_coupling_mean_no_coherent(const BaseScatterer &scatterer) const
{
    return get_coupling_point_no_coherent(scatterer);
}

double Detector::get_coupling_point_coherent(const BaseScatterer &scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    auto [horizontal_projection, vertical_projection] = this->get_projected_farfields(theta_field, phi_field);

    this->apply_scalar_field(horizontal_projection, vertical_projection);

    double
        coupling_theta = get_norm1_squared(horizontal_projection),
        coupling_phi = get_norm1_squared(vertical_projection);

    this->apply_polarization_filter(
        coupling_theta,
        coupling_phi,
        this->polarization_filter
    );

    return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->fibonacci_mesh.dOmega;
}

double Detector::get_coupling_point_no_coherent(const BaseScatterer &scatterer) const
{
    auto [theta_field, phi_field] = scatterer.compute_unstructured_farfields(this->fibonacci_mesh, 1.0);

    double
        coupling_theta = this->get_norm2_squared(theta_field),
        coupling_phi = this->get_norm2_squared(phi_field);

    this->apply_polarization_filter(
        coupling_theta,
        coupling_phi,
        this->polarization_filter
    );

    return 0.5 * EPSILON0 * C * (coupling_theta + coupling_phi) * this->fibonacci_mesh.dOmega;
}
