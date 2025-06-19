#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "fibonacci.h"

void register_fibonacci(pybind11::module_& module) {
    pybind11::class_<FibonacciMesh>(module, "FIBONACCIMESH")
        .def(pybind11::init<size_t, double, double, double, double, double>(),
            pybind11::arg("sampling"),
            pybind11::arg("max_angle"),
            pybind11::arg("min_angle"),
            pybind11::arg("phi_offset"),
            pybind11::arg("gamma_offset"),
            pybind11::arg("rotation_angle"),
            R"pbdoc(
                Initialize a Fibonacci mesh with specified parameters.

                Parameters
                ----------
                sampling : int
                    Number of sampling points in the mesh.
                min_angle : float
                    Minimum angle for the mesh points (in radians).
                max_angle : float
                    Maximum angle for the mesh points (in radians).
                phi_offset : float
                    Offset for the azimuthal angle (phi) in radians.
                rotation_angle : float
                    Rotation angle for the mesh in radians.
                gamma_offset : float
                    Offset for the polar angle (gamma) in radians.
            )pbdoc"
        )
        .def_readonly("_cpp_vertical_base",
            &FibonacciMesh::vertical_vector_field,
            R"pbdoc(
                Vertical base vector field of the Fibonacci mesh.
                This property provides access to the vertical vector field used in the mesh.
            )pbdoc"
        )
        .def_readonly("_cpp_horizontal_base",
            &FibonacciMesh::horizontal_vector_field,
            R"pbdoc(
                Horizontal base vector field of the Fibonacci mesh.
                This property provides access to the horizontal vector field used in the mesh.
            )pbdoc"
        )
        .def_readonly("cartesian",
            &FibonacciMesh::cartesian,
            R"pbdoc(
                Cartesian coordinates of the Fibonacci mesh points.
                This property provides access to the x, y, and z coordinates of the mesh points.
            )pbdoc"
        )
        .def_readonly("spherical",
            &FibonacciMesh::spherical,
            R"pbdoc(
                Spherical coordinates of the Fibonacci mesh points.
                This property provides access to the radial distance, azimuthal angle, and polar angle of the mesh points.
            )pbdoc"
        )
        .def_readonly("base_cartesian",
            &FibonacciMesh::base_cartesian,
            R"pbdoc(
                Base Cartesian coordinates of the Fibonacci mesh points.
                This property provides access to the base x, y, and z coordinates of the mesh points.
            )pbdoc"
        )
        // Differential solid angle and total solid angle
        .def_readwrite("_cpp_d_omega",
            &FibonacciMesh::dOmega,
            "Differential solid angle covered by each point."
        )
        .def_readwrite("_cpp_omega",
            &FibonacciMesh::Omega,
            "Total solid angle covered by the mesh."
        )
        .def_readonly("_cpp_phi_offset",
            &FibonacciMesh::phi_offset,
            "Azimuthal angle offset (in radians) for the mesh."
        )
        .def_readonly("_cpp_gamma_offset",
            &FibonacciMesh::gamma_offset,
            "Polar angle offset (in radians) for the mesh."
        )
        .def_readonly("_cpp_rotation",
            &FibonacciMesh::rotation,
            "Rotation angle (in radians) for the mesh."
        )
        // Methods for rotation and field computation
        .def("rotate_around_axis",
            &FibonacciMesh::rotate_around_axis,
            pybind11::arg("angle"),
            "Rotates the mesh around a specified axis by a given angle."
        )
        .def("compute_vector_field",
            &FibonacciMesh::compute_vector_field,
            "Computes the vector field based on the mesh configuration."
        )
        .def("compute_projections",
            &FibonacciMesh::compute_projections,
            "Computes projections of the vector field."
        )
        // Properties for vector field projections
        .def_readonly("perpendicular",
            &FibonacciMesh::perpendicular_vector,
            "Perpendicular vector field component."
        )
        .def_readonly("parallel",
            &FibonacciMesh::parallel_vector,
            "Parallel vector field component."
        )
        // Properties for vector field projections along H and V directions
        .def_property_readonly("horizontal_to_parallel",
            [](const FibonacciMesh& self) {
                return pybind11::array_t<double>(self.horizontal_parallel_projection.size(), self.horizontal_parallel_projection.data());
            },
            R"pbdoc(
                Returns the projection of the horizontal base vector onto the parallel component.

                Returns
                -------
                numpy.ndarray
                    The horizontal projection onto the parallel vector field.
            )pbdoc"
        )
        .def_property_readonly("horizontal_to_perpendicular",
            [](const FibonacciMesh& self) {
                return pybind11::array_t<double>(self.horizontal_perpendicular_projection.size(), self.horizontal_perpendicular_projection.data());
            },
            R"pbdoc(
                Returns the projection of the horizontal base vector onto the perpendicular component.

                Returns
                -------
                numpy.ndarray
                    The horizontal projection onto the perpendicular vector field.
            )pbdoc"
        )
        .def_property_readonly("vertical_to_parallel",
            [](const FibonacciMesh& self) {
                return pybind11::array_t<double>(self.vertical_parallel_projection.size(), self.vertical_parallel_projection.data());
            },
            R"pbdoc(
                Returns the projection of the vertical base vector onto the parallel component.

                Returns
                -------
                numpy.ndarray
                    The vertical projection onto the parallel vector field.
            )pbdoc"
        )
        .def_property_readonly("vertical_to_perpendicular",
            [](const FibonacciMesh& self) {
                return pybind11::array_t<double>(self.vertical_perpendicular_projection.size(), self.vertical_perpendicular_projection.data());
            },
            R"pbdoc(
                Returns the projection of the vertical base vector onto the perpendicular component.

                Returns
                -------
                numpy.ndarray
                    The vertical projection onto the perpendicular vector field.
            )pbdoc"
        )
    ;
}

