#include "fibonacci/fibonacci.h"

#define PI (double)3.14159265358979323846264338

void FibonacciMesh::rotate_around_center() {
    if (gamma_offset != 0.0) {
        cartesian_coordinates.rotate_about_axis('x', gamma_offset);
        vertical_vector_field.rotate_about_axis('x', gamma_offset);
        horizontal_vector_field.rotate_about_axis('x', gamma_offset);
    }

    if (phi_offset != 0.0) {
        cartesian_coordinates.rotate_about_axis('y', phi_offset);
        vertical_vector_field.rotate_about_axis('y', phi_offset);
        horizontal_vector_field.rotate_about_axis('y', phi_offset);
    }
}

void FibonacciMesh::compute_vector_field(){
    parallel_vector = VectorField(sampling);
    perpendicular_vector = VectorField(sampling);

    for (size_t i = 0; i < sampling; i++){
        this->perpendicular_vector.at(i,0) = -cartesian_coordinates.y[i] ;
        this->perpendicular_vector.at(i,1) = cartesian_coordinates.x[i] ;
        this->perpendicular_vector.at(i,2) = 0.0;

        this->parallel_vector.at(i,0) = cartesian_coordinates.x[i] * cartesian_coordinates.z[i] ;
        this->parallel_vector.at(i,1) = cartesian_coordinates.y[i] * cartesian_coordinates.z[i] ;
        this->parallel_vector.at(i,2) = -( pow( cartesian_coordinates.x[i], 2 ) + pow( cartesian_coordinates.y[i], 2 ) );
    }

    parallel_vector.normalize();
    perpendicular_vector.normalize();
}

void FibonacciMesh::compute_projections(){
    horizontal_parallel_projection = parallel_vector.get_scalar_product(horizontal_vector_field);

    vertical_parallel_projection = parallel_vector.get_scalar_product(vertical_vector_field);

    horizontal_perpendicular_projection = perpendicular_vector.get_scalar_product(horizontal_vector_field);

    vertical_perpendicular_projection = perpendicular_vector.get_scalar_product(vertical_vector_field);
}

void FibonacciMesh::compute_mesh(){
    this->compute_properties();

    double golden_angle = PI * (3. - sqrt(5.));  // golden angle = 2.39996322972865332

    for (size_t i = 0; i < true_number_of_sample; i++){
        double
            z = 1 - (2. * i) / (true_number_of_sample - 1),
            theta = golden_angle * i,
            radius = sqrt(1 - z * z),
            angle = atan(radius / z);

        if (angle < this->min_angle)
            continue;
        if (cartesian_coordinates.z.size() >= this->sampling)
            break;
        cartesian_coordinates.z.push_back(z);
        cartesian_coordinates.x.push_back(cos(theta) * radius);
        cartesian_coordinates.y.push_back(sin(theta) * radius);
    }
}


void FibonacciMesh::compute_properties(){
    double
        solid_angle = 2 * PI * abs(cos(min_angle) - cos(max_angle)),
        ratio = ( 4. * PI / solid_angle );

    this->Omega = solid_angle;
    this->dOmega = Omega / this->sampling;
    this->true_number_of_sample = (size_t) ( this->sampling * ratio );
}

std::vector<double> FibonacciMesh::get_principal_axis() const {
    std::vector<double> principal_axis = {
        cartesian_coordinates.x[0],
        cartesian_coordinates.y[0],
        cartesian_coordinates.z[0]
    };

    return principal_axis;
}

void FibonacciMesh::rotate_around_axis(double angle) {
    std::vector<double> rotation_axis = this->get_principal_axis();

    std::vector<std::vector<double>> rotation_matrix = get_rotation_matrix(rotation_axis, angle);

    for (size_t i = 0; i < cartesian_coordinates.x.size(); ++i){
        double
            x = cartesian_coordinates.x[i],
            y = cartesian_coordinates.y[i],
            z = cartesian_coordinates.z[i];

        cartesian_coordinates.x[i] = rotation_matrix[0][0] * x + rotation_matrix[0][1] * y + rotation_matrix[0][2] * z;
        cartesian_coordinates.y[i] = rotation_matrix[1][0] * x + rotation_matrix[1][1] * y + rotation_matrix[1][2] * z;
        cartesian_coordinates.z[i] = rotation_matrix[2][0] * x + rotation_matrix[2][1] * y + rotation_matrix[2][2] * z;
    }

    this->spherical_coordinates = this->cartesian_coordinates.to_spherical();
}

// --
