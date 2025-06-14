#pragma once

#define DEFINE_PY_GETTER(dtype, py_name, name) \
    py::array_t<dtype> get_##py_name() const {return _vector_to_numpy(name, {name.size()}); };


#include <vector>
#include <cmath>
#include "utils/coordinates.cpp"
#include "utils/numpy_interface.h"
#include "utils/base_mesh.h"

namespace py = pybind11;

class FibonacciMesh : public BaseMesh {
    public:
        size_t true_number_of_sample = 0;
        double max_angle = 0.0;
        double min_angle = 0.0;
        double phi_offset = 0.0;
        double gamma_offset = 0.0;
        double dOmega = 0.0;
        double Omega = 0.0;

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

        FibonacciMesh(int sampling, double max_angle, double min_angle, double phi_offset, double gamma_offset, double rotation, double radius = 1.0):
            BaseMesh(sampling, radius), max_angle(max_angle), min_angle(min_angle), phi_offset(phi_offset), gamma_offset(gamma_offset) {

            this->compute_mesh();
            base_cartesian_coordinates = cartesian_coordinates;

            this->rotate_around_center();
            this->rotate_around_axis(rotation);

            this->compute_vector_field();
            this->compute_projections();
        }

        void rotate_around_center();
        void compute_vector_field();
        void compute_projections();
        void rotate_around_axis(double angle);
        void compute_mesh();
        void compute_properties();

        std::vector<double> get_principal_axis() const;

        py::array_t<double> get_parallel_vector() const {
            return _vector_to_numpy(parallel_vector.data, {parallel_vector.sampling, 3});
        };
        py::array_t<double> get_perpendicular_vector() const {
            return _vector_to_numpy(perpendicular_vector.data, {perpendicular_vector.sampling, 3});
        };

        std::vector<std::vector<double>> get_rotation_matrix(std::vector<double> rotation_axis, double rotation_angle) const
        {
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

        DEFINE_PY_GETTER(double, horizontal_parallel_projection, horizontal_parallel_projection)
        DEFINE_PY_GETTER(double, vertical_parallel_projection, vertical_parallel_projection)
        DEFINE_PY_GETTER(double, horizontal_perpendicular_projection, horizontal_perpendicular_projection)
        DEFINE_PY_GETTER(double, vertical_perpendicular_projection, vertical_perpendicular_projection)

        DEFINE_PY_GETTER(double, x_py, cartesian_coordinates.x)
        DEFINE_PY_GETTER(double, y_py, cartesian_coordinates.y)
        DEFINE_PY_GETTER(double, z_py, cartesian_coordinates.z)

        DEFINE_PY_GETTER(double, base_x_py, base_cartesian_coordinates.x)
        DEFINE_PY_GETTER(double, base_y_py, base_cartesian_coordinates.y)
        DEFINE_PY_GETTER(double, base_z_py, base_cartesian_coordinates.z)

        DEFINE_PY_GETTER(double, r_py, spherical_coordinates.r)
        DEFINE_PY_GETTER(double, phi_py, spherical_coordinates.phi)
        DEFINE_PY_GETTER(double, theta_py, spherical_coordinates.theta)
};
