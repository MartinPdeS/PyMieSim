#include <pybind11/pybind11.h>
#include "fibonnaci_mesh.cpp"

PYBIND11_MODULE(Fibonacci, module) {
    module.doc() = "LGeneralized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";

    pybind11::class_<FibonacciMesh>(module, "FIBONACCI")
        .def(
            pybind11::init<int, double, double, double, double>(),
            pybind11::arg("sampling"),
            pybind11::arg("max_angle"),
            pybind11::arg("phi_offset"),
            pybind11::arg("rotation_angle"),
            pybind11::arg("gamma_offset")
        )

        .def_property("x", &FibonacciMesh::get_x_py, &FibonacciMesh::set_x_py)
        .def_property("y", &FibonacciMesh::get_y_py, &FibonacciMesh::set_y_py)
        .def_property("z", &FibonacciMesh::get_z_py, &FibonacciMesh::set_z_py)

        .def_property_readonly("r", &FibonacciMesh::get_r_py)
        .def_property_readonly("phi", &FibonacciMesh::get_phi_py)
        .def_property_readonly("theta", &FibonacciMesh::get_theta_py)

        .def_readwrite("d_omega", &FibonacciMesh::dOmega)
        .def_readwrite("omega", &FibonacciMesh::Omega)

        .def("rotate_around_axis", &FibonacciMesh::rotate_around_axis)
        .def("compute_vector_field", &FibonacciMesh::compute_vector_field)
        .def("compute_projections", &FibonacciMesh::compute_projections)

        .def_property_readonly("parallel_vector", &FibonacciMesh::get_parallel_vector)
        .def_property_readonly("perpendicular_vector", &FibonacciMesh::get_perpendicular_vector)

        .def_property_readonly("H_para", &FibonacciMesh::get_horizontal_parallel_projection)
        .def_property_readonly("H_perp", &FibonacciMesh::get_horizontal_perpendicular_projection)
        .def_property_readonly("V_para", &FibonacciMesh::get_vertical_parallel_projection)
        .def_property_readonly("V_perp", &FibonacciMesh::get_vertical_perpendicular_projection)
        ;
  }







// -
