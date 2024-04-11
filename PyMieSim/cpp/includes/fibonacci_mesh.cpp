#pragma once

#include "coordinates.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"
#include "fibonacci_mesh.h"




void FibonacciMesh::rotate_around_center() {
    if (gamma_offset != 0.0) {
        cartesian_coordinates.mx_rot_x(gamma_offset);
        vertical_vector_field.mx_rot_x(gamma_offset);
        horizontal_vector_field.mx_rot_x(gamma_offset);
    }

    if (phi_offset != 0.0) {
        cartesian_coordinates.mx_rot_y(phi_offset);
        vertical_vector_field.mx_rot_y(phi_offset);
        horizontal_vector_field.mx_rot_y(phi_offset);
    }
}

void FibonacciMesh::compute_vector_field(){
    parallel_vector = VectorField(sampling);
    perpendicular_vector = VectorField(sampling);

    for (size_t i = 0; i < sampling; i++){
        this->perpendicular_vector(i,0) = -cartesian_coordinates.Y[i] ;
        this->perpendicular_vector(i,1) = cartesian_coordinates.X[i] ;
        this->perpendicular_vector(i,2) = 0.0;

        this->parallel_vector(i,0) = cartesian_coordinates.X[i] * cartesian_coordinates.Z[i] ;
        this->parallel_vector(i,1) = cartesian_coordinates.Y[i] * cartesian_coordinates.Z[i] ;
        this->parallel_vector(i,2) = -( pow( cartesian_coordinates.X[i], 2 ) + pow( cartesian_coordinates.Y[i], 2 ) );
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

double FibonacciMesh::NA2Angle(const double &NA) const {
    if (NA <= 1.0)
    return asin(NA);

    if (NA >= 1.0)
    return asin(NA-1.0) + PI/2.0;

    return 1.0;
}


void FibonacciMesh::compute_mesh(){
    this->compute_properties();

    double golden_angle = PI * (3. - sqrt(5.));  // golden angle = 2.39996322972865332

    for (size_t i = 0; i < true_number_of_sample; i++){
        cartesian_coordinates.Z[i] = 1 - (i / (double)( true_number_of_sample - 1) ) * 2 ;

        double
            theta = golden_angle * i,
            radius = sqrt(1 - cartesian_coordinates.Z[i] * cartesian_coordinates.Z[i]);

        cartesian_coordinates.X[i] = cos(theta) * radius ;
        cartesian_coordinates.Y[i] = sin(theta) * radius ;

        if (i == (size_t) sampling - 1)
            break;
    }
}

void FibonacciMesh::compute_properties(){
    double
        solid_angle = abs( 2. * PI * ( cos(max_angle) - 1. ) ),   //cos(0) =1
        ratio = ( 4. * PI / solid_angle );

    size_t _true_number_of_sample = (size_t) ( sampling * ratio );

    dOmega = 4. * PI / _true_number_of_sample;
    Omega = dOmega * sampling ;
    true_number_of_sample = _true_number_of_sample;
}

std::vector<double> FibonacciMesh::get_principal_axis() const {
    std::vector<double> principal_axis = {
        cartesian_coordinates.X[0],
        cartesian_coordinates.Y[0],
        cartesian_coordinates.Z[0]
    };

    return principal_axis;
}

void FibonacciMesh::rotate_around_axis(double angle) {
    std::vector<double> rotation_axis = this->get_principal_axis();

    std::vector<std::vector<double>> rotation_matrix = get_rotation_matrix(rotation_axis, angle);

    for (size_t i = 0; i < cartesian_coordinates.X.size(); ++i){
        double
            x = cartesian_coordinates.X[i],
            y = cartesian_coordinates.Y[i],
            z = cartesian_coordinates.Z[i];

        cartesian_coordinates.X[i] = rotation_matrix[0][0] * x + rotation_matrix[0][1] * y + rotation_matrix[0][2] * z;
        cartesian_coordinates.Y[i] = rotation_matrix[1][0] * x + rotation_matrix[1][1] * y + rotation_matrix[1][2] * z;
        cartesian_coordinates.Z[i] = rotation_matrix[2][0] * x + rotation_matrix[2][1] * y + rotation_matrix[2][2] * z;
    }

    this->spherical_coordinates = this->cartesian_coordinates.cartesian_to_spherical();
}


class FullSteradian
{
    public:
        size_t sampling;

        Spherical spherical_coordinates;

        double dTheta, dPhi;

        FullSteradian(const size_t sampling) : sampling(sampling), spherical_coordinates(sampling)
        {
            dTheta = 2.0 * PI / (sampling-1);
            dPhi   = 1.0 * PI / (sampling-1);

            for (size_t p=0; p<sampling; p++)
                spherical_coordinates.Phi[p]   = p * dPhi   - PI/2.0;

            for (size_t t=0; t<sampling; t++)
                spherical_coordinates.Theta[t] = t * dTheta - PI/1.0;
        }

        template<typename T>
        T get_integral(std::vector<T>& Vector)
        {
            T integral = 0;

            for (size_t p=0; p<sampling; p++)
                for (size_t t=0; t<sampling; t++)
                    integral += Vector[p*sampling + t] * sin(spherical_coordinates.Phi[p] + PI/2.0) * dPhi * dTheta;

            return integral;
        }


        template<typename T>
        T get_cos_integral(const std::vector<T>& vector) const
        {
            T integral = 0;

            for (size_t p=0; p<sampling; p++)
                for (size_t t=0; t<sampling; t++)
                    integral += vector[p*sampling + t] * cos(spherical_coordinates.Phi[p] + PI/2.0) * sin(spherical_coordinates.Phi[p] + PI/2.0) * dPhi * dTheta;

            return integral;
        }

        double get_integral() const
        {
            double integral = 0;

            for (auto phi : spherical_coordinates.Phi)
                for (auto theta : spherical_coordinates.Theta)
                    integral += sin(phi+PI/2.0) * dPhi * dTheta;

            return integral;
        }

};

// --
