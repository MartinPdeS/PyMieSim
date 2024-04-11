#ifndef FIBONACCI_MESH_H
#define FIBONACCI_MESH_H

#include <vector>
#include <cmath>
#include "coordinates.cpp"
#include "utils.cpp"
#include "numpy_interface.cpp"

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

    FibonacciMesh(int sampling, double max_angle, double phi_offset,
        double gamma_offset, double rotation_angle):
        sampling(sampling), max_angle(max_angle),
        phi_offset(phi_offset), gamma_offset(gamma_offset) {

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

    py::array_t<double> get_parallel_vector() const {return parallel_vector.get_numpy();};
    py::array_t<double> get_perpendicular_vector() const {return perpendicular_vector.get_numpy();};

    py::array_t<double> get_horizontal_parallel_projection() const {return vector_to_numpy_copy(horizontal_parallel_projection);};
    py::array_t<double> get_vertical_parallel_projection() const {return vector_to_numpy_copy(vertical_parallel_projection);};
    py::array_t<double> get_horizontal_perpendicular_projection() const {return vector_to_numpy_copy(horizontal_perpendicular_projection);};
    py::array_t<double> get_vertical_perpendicular_projection() const {return vector_to_numpy_copy(vertical_perpendicular_projection);};

    py::array_t<double> get_x_py() const {return cartesian_coordinates.get_x_py();};
    py::array_t<double> get_y_py() const {return cartesian_coordinates.get_y_py();};
    py::array_t<double> get_z_py() const {return cartesian_coordinates.get_z_py();};

    py::array_t<double> get_base_x_py() const {return base_cartesian_coordinates.get_x_py();};
    py::array_t<double> get_base_y_py() const {return base_cartesian_coordinates.get_y_py();};
    py::array_t<double> get_base_z_py() const {return base_cartesian_coordinates.get_z_py();};

    void set_x_py(const std::vector<double> &value) {cartesian_coordinates.set_x_py(value);};
    void set_y_py(const std::vector<double> &value) {cartesian_coordinates.set_y_py(value);};
    void set_z_py(const std::vector<double> &value) {cartesian_coordinates.set_z_py(value);};

    py::array_t<double> get_r_py() const {return spherical_coordinates.get_r_py();};
    py::array_t<double> get_phi_py() const {return spherical_coordinates.get_phi_py();};
    py::array_t<double> get_theta_py() const { return spherical_coordinates.get_theta_py();};

};

#endif // FIBONACCI_MESH_H
