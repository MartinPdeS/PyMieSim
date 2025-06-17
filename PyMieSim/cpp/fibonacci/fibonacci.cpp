#include "fibonacci/fibonacci.h"


// ------------------ Constructors ------------------
FibonacciMesh::FibonacciMesh(size_t sampling, double max_angle, double min_angle, double phi_offset, double gamma_offset, double rotation, double radius):
    sampling(sampling), max_angle(max_angle), min_angle(min_angle), phi_offset(phi_offset), gamma_offset(gamma_offset), radius(radius), cartesian(sampling), spherical(sampling)
{
    this->compute_properties();
    this->compute_mesh();
    this->base_cartesian = cartesian;

    this->rotate_around_center();
    this->rotate_around_axis(rotation);

    this->compute_vector_field();
    this->compute_projections();
}

// ------------------ Methods ------------------
void FibonacciMesh::rotate_around_center() {
    if (gamma_offset != 0.0) {
        this->cartesian.rotate_about_axis('x', this->gamma_offset);
        this->vertical_vector_field.rotate_about_axis('x', this->gamma_offset);
        this->horizontal_vector_field.rotate_about_axis('x', this->gamma_offset);
    }

    if (phi_offset != 0.0) {
        this->cartesian.rotate_about_axis('y', this->phi_offset);
        this->vertical_vector_field.rotate_about_axis('y', this->phi_offset);
        this->horizontal_vector_field.rotate_about_axis('y', this->phi_offset);
    }
}

std::vector<std::vector<double>>
FibonacciMesh::get_rotation_matrix(std::vector<double> rotation_axis, double rotation_angle) const {
    double norm_rotation_axis = sqrt(pow(rotation_axis[0], 2) + pow(rotation_axis[1], 2) + pow(rotation_axis[2], 2));

    for (double &x: rotation_axis)
        x /= norm_rotation_axis;

    double
        a = cos(rotation_angle / 2.0),
        b = -1 * sin(rotation_angle / 2.0) * rotation_axis[0],
        c = -1 * sin(rotation_angle / 2.0) * rotation_axis[1],
        d = -1 * sin(rotation_angle / 2.0) * rotation_axis[2];

    std::vector<std::vector<double>> matrix = {
        {a * a + b * b - c * c - d * d, 2 * (b * c + a * d), 2 * (b * d - a * c)},
        {2 * (b * c - a * d), a * a + c * c - b * b - d * d, 2 * (c * d + a * b)},
        {2 * (b * d + a * c), 2 * (c * d - a * b), a * a + d * d - b * b - c * c}
    };

    return matrix;
}

void FibonacciMesh::compute_vector_field() {
    parallel_vector = VectorField(sampling);
    perpendicular_vector = VectorField(sampling);

    for (size_t i = 0; i < sampling; i++){
        this->perpendicular_vector.at(i,0) = -cartesian.y[i] ;
        this->perpendicular_vector.at(i,1) = cartesian.x[i] ;
        this->perpendicular_vector.at(i,2) = 0.0;

        this->parallel_vector.at(i,0) = this->cartesian.x[i] * this->cartesian.z[i] ;
        this->parallel_vector.at(i,1) = this->cartesian.y[i] * this->cartesian.z[i] ;
        this->parallel_vector.at(i,2) = -( pow( this->cartesian.x[i], 2 ) + pow( this->cartesian.y[i], 2 ) );
    }

    this->parallel_vector.normalize();
    this->perpendicular_vector.normalize();
}

void FibonacciMesh::compute_projections() {
    this->horizontal_parallel_projection = this->parallel_vector.get_scalar_product(this->horizontal_vector_field);

    this->vertical_parallel_projection = this->parallel_vector.get_scalar_product(this->vertical_vector_field);

    this->horizontal_perpendicular_projection = this->perpendicular_vector.get_scalar_product(this->horizontal_vector_field);

    this->vertical_perpendicular_projection = this->perpendicular_vector.get_scalar_product(this->vertical_vector_field);
}

void FibonacciMesh::compute_mesh() {
    double golden_angle = PI * (3. - sqrt(5.));  // golden angle = 2.39996322972865332

    for (size_t i = 0; i < this->true_number_of_sample; i++){
        double
            z = 1 - (2. * i) / (true_number_of_sample - 1),
            theta = golden_angle * i,
            radius = sqrt(1 - z * z),
            angle = atan(radius / z);

        if (angle < this->min_angle)
            continue;

        if (this->cartesian.z.size() >= this->sampling)
            break;

        this->cartesian.z.push_back(z);
        this->cartesian.x.push_back(cos(theta) * radius);
        this->cartesian.y.push_back(sin(theta) * radius);
    }
}

void FibonacciMesh::compute_properties(){
    double
        diff = cos(this->min_angle) - cos(this->max_angle),
        solid_angle = 2. * PI * std::abs(diff),
        ratio = ( 4. * PI / solid_angle );

    this->Omega = solid_angle;
    this->dOmega = Omega / this->sampling;
    this->true_number_of_sample = (size_t) ( this->sampling * ratio );
}

std::vector<double> FibonacciMesh::get_principal_axis() const {
    std::vector<double> principal_axis = {
        this->cartesian.x[0],
        this->cartesian.y[0],
        this->cartesian.z[0]
    };

    return principal_axis;
}

void FibonacciMesh::rotate_around_axis(double angle) {
    std::vector<double> rotation_axis = this->get_principal_axis();

    std::vector<std::vector<double>> rotation_matrix = this->get_rotation_matrix(rotation_axis, angle);

    for (size_t i = 0; i < cartesian.x.size(); ++i){
        double
            x = this->cartesian.x[i],
            y = this->cartesian.y[i],
            z = this->cartesian.z[i];

        this->cartesian.x[i] = rotation_matrix[0][0] * x + rotation_matrix[0][1] * y + rotation_matrix[0][2] * z;
        this->cartesian.y[i] = rotation_matrix[1][0] * x + rotation_matrix[1][1] * y + rotation_matrix[1][2] * z;
        this->cartesian.z[i] = rotation_matrix[2][0] * x + rotation_matrix[2][1] * y + rotation_matrix[2][2] * z;
    }

    this->spherical = this->cartesian.to_spherical();
}

// --
