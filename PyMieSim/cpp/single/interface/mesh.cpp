#include <pybind11/pybind11.h>
#include "single/includes/fibonacci_mesh.cpp"

namespace py = pybind11;

PYBIND11_MODULE(Fibonacci, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) C++ binding module for light scattering from a spherical scatterer";

    py::class_<FibonacciMesh>(module, "FibonacciMesh")
        .def(py::init<int, double, double, double, double, double>(),
             py::arg("sampling"),
             py::arg("min_angle") = 0.,
             py::arg("max_angle"),
             py::arg("phi_offset"),
             py::arg("rotation_angle"),
             py::arg("gamma_offset"),
             "Initializes a Fibonacci mesh with specified parameters.")

        // Properties for coordinates
        .def_property_readonly("x", &FibonacciMesh::get_x_py, "X coordinates of points on the Fibonacci mesh.")
        .def_property_readonly("y", &FibonacciMesh::get_y_py, "Y coordinates of points on the Fibonacci mesh.")
        .def_property_readonly("z", &FibonacciMesh::get_z_py, "Z coordinates of points on the Fibonacci mesh.")

        .def_property_readonly("base_x", &FibonacciMesh::get_base_x_py, "X coordinates of points on the Fibonacci mesh.")
        .def_property_readonly("base_y", &FibonacciMesh::get_base_y_py, "Y coordinates of points on the Fibonacci mesh.")
        .def_property_readonly("base_z", &FibonacciMesh::get_base_z_py, "Z coordinates of points on the Fibonacci mesh.")

        // Readonly properties for spherical coordinates and differential solid angles
        .def_property_readonly("r", &FibonacciMesh::get_r_py, "Radial distance of points on the Fibonacci mesh.")
        .def_property_readonly("phi", &FibonacciMesh::get_phi_py, "Azimuthal angle (phi) of points on the Fibonacci mesh.")
        .def_property_readonly("theta", &FibonacciMesh::get_theta_py, "Polar angle (theta) of points on the Fibonacci mesh.")

        // Differential solid angle and total solid angle
        .def_readwrite("d_omega", &FibonacciMesh::dOmega, "Differential solid angle covered by each point.")
        .def_readwrite("omega", &FibonacciMesh::Omega, "Total solid angle covered by the mesh.")

        // Methods for rotation and field computation
        .def("rotate_around_axis", &FibonacciMesh::rotate_around_axis, py::arg("angle"), "Rotates the mesh around a specified axis by a given angle.")
        .def("compute_vector_field", &FibonacciMesh::compute_vector_field, "Computes the vector field based on the mesh configuration.")
        .def("compute_projections", &FibonacciMesh::compute_projections, "Computes projections of the vector field.")

        // Properties for vector field projections
        .def_property_readonly("parallel_vector", &FibonacciMesh::get_parallel_vector, "Parallel component of the vector field.")
        .def_property_readonly("perpendicular_vector", &FibonacciMesh::get_perpendicular_vector, "Perpendicular component of the vector field.")

        // Properties for vector field projections along H and V directions
        .def_property_readonly("H_para", &FibonacciMesh::get_horizontal_parallel_projection, "Horizontal parallel projection of the vector field.")
        .def_property_readonly("H_perp", &FibonacciMesh::get_horizontal_perpendicular_projection, "Horizontal perpendicular projection of the vector field.")
        .def_property_readonly("V_para", &FibonacciMesh::get_vertical_parallel_projection, "Vertical parallel projection of the vector field.")
        .def_property_readonly("V_perp", &FibonacciMesh::get_vertical_perpendicular_projection, "Vertical perpendicular projection of the vector field.");
}
