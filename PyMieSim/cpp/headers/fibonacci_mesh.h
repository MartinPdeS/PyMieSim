#ifndef FIBONACCI_MESH_H
#define FIBONACCI_MESH_H

#include <vector>
#include <cmath>
#include "coordinates.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"

namespace py = pybind11;

class FibonacciMesh {
    public:
        size_t sampling = 0;
        size_t true_number_of_sample = 0;
        double max_angle = 0.0;
        double phi_offset = 0.0;
        double gamma_offset = 0.0;
        double dOmega = 0.0;
        double Omega = 0.0;

        Spherical spherical_coordinates;
        Cartesian cartesian_coordinates;
        Cartesian base_cartesian_coordinates;

        VectorField perpendicular_vector;
        VectorField parallel_vector;
        VectorField vertical_vector_field{{1, 0, 0}};
        VectorField horizontal_vector_field{{0, 1, 0}};

        std::vector<double> horizontal_parallel_projection;
        std::vector<double> vertical_parallel_projection;
        std::vector<double> horizontal_perpendicular_projection;
        std::vector<double> vertical_perpendicular_projection;

        FibonacciMesh() = default;

        FibonacciMesh(int sampling, double max_angle, double phi_offset, double gamma_offset, double rotation_angle):
            sampling(sampling), max_angle(max_angle), phi_offset(phi_offset), gamma_offset(gamma_offset) {

            cartesian_coordinates = Cartesian(sampling);
            compute_mesh();
            base_cartesian_coordinates = cartesian_coordinates;

            rotate_around_center();
            rotate_around_axis(rotation_angle);

            compute_vector_field();
            compute_projections();
        }

        void rotate_around_center();
        void compute_vector_field();
        void compute_projections();
        void rotate_around_axis(double angle);
        void compute_mesh();
        void compute_properties();

        double NA2Angle(const double &NA) const;

        std::vector<double> get_principal_axis() const;

        py::array_t<double> get_parallel_vector() const {return vector_to_numpy_copy(parallel_vector.data, {parallel_vector.sampling, 3});};
        py::array_t<double> get_perpendicular_vector() const {return vector_to_numpy_copy(perpendicular_vector.data, {perpendicular_vector.sampling, 3});};

        py::array_t<double> get_horizontal_parallel_projection() const {return vector_to_numpy_copy(horizontal_parallel_projection);};
        py::array_t<double> get_vertical_parallel_projection() const {return vector_to_numpy_copy(vertical_parallel_projection);};
        py::array_t<double> get_horizontal_perpendicular_projection() const {return vector_to_numpy_copy(horizontal_perpendicular_projection);};
        py::array_t<double> get_vertical_perpendicular_projection() const {return vector_to_numpy_copy(vertical_perpendicular_projection);};

        py::array_t<double> get_x_py() const {return vector_to_numpy_copy(cartesian_coordinates.x);};
        py::array_t<double> get_y_py() const {return vector_to_numpy_copy(cartesian_coordinates.y);};
        py::array_t<double> get_z_py() const {return vector_to_numpy_copy(cartesian_coordinates.z);};

        py::array_t<double> get_base_x_py() const {return vector_to_numpy_copy(base_cartesian_coordinates.x);};
        py::array_t<double> get_base_y_py() const {return vector_to_numpy_copy(base_cartesian_coordinates.y);};
        py::array_t<double> get_base_z_py() const {return vector_to_numpy_copy(base_cartesian_coordinates.z);};

        py::array_t<double> get_r_py() const {return vector_to_numpy_copy(spherical_coordinates.r);};
        py::array_t<double> get_phi_py() const {return vector_to_numpy_copy(spherical_coordinates.phi);};
        py::array_t<double> get_theta_py() const { return vector_to_numpy_copy(spherical_coordinates.theta);};

};

#endif // FIBONACCI_MESH_H
