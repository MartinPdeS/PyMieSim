#include <pybind11/pybind11.h>
#include "detector/detector.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_detector, module) {
    module.doc() = R"pbdoc(
        Detector binding for PyMieSim.

        Provides a `BindedDetector` class to compute coupling to scatterers
        and expose the detected field distribution.
    )pbdoc";

    py::class_<FibonacciMesh>(module, "FIBONACCIMESH")
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
        .def_property_readonly("V_perp", &FibonacciMesh::get_vertical_perpendicular_projection, "Vertical perpendicular projection of the vector field.")
        ;


    py::class_<ModeField>(module, "MODEFIELD")
        .def("_cpp_get_unstructured",
            &ModeField::get_unstructured,
            py::arg("x_coords"),
            py::arg("y_coords"),
            py::return_value_policy::move,
            R"pbdoc(
                Generate a structured scalar field as a numpy array.

                Parameters
                ----------
                sampling : int
                    The sampling rate for the scalar field. Default is 100.

                Returns
                -------
                numpy.ndarray
                    A 2D array representing the structured scalar field.
            )pbdoc")

    ;

    py::class_<Detector>(module, "DETECTOR", R"pbdoc(
        Detector for Lorenz-Mie scattering simulations.

        This class wraps a C++ Detector that samples the scattered field
        according to numerical aperture, orientation offsets, polarization
        filtering, and optional coherent or mean coupling strategies.

        Parameters
        ----------
        mode_number : str
            Identifier for the detector mode (e.g. '0', '1', ...).
        sampling : int
            Number of sample points on the detector's angular mesh.
        NA : float
            Numerical aperture of the detector; controls the angular acceptance.
        cache_NA : float
            Cached NA threshold for fast lookup (set equal to `NA` to disable).
        phi_offset : float
            Azimuthal angle offset (radians) for detector orientation.
        gamma_offset : float
            Polar angle offset (radians) for detector orientation.
        polarization_filter : float
            Degree or type of polarization filtering applied.
        rotation : float
            Rotation angle (radians) of the detector's field-of-view.
        coherent : bool
            Whether to compute coherent coupling (true) or incoherent (false).
        mean_coupling : bool
            If true, uses mean coupling; if false, uses point coupling.
        medium_refractive_index : float, optional
            Refractive index of the surrounding medium. Default is 1.0.

        Attributes
        ----------
        scalar_field : numpy.ndarray
            The detected scalar field (intensity) sampled on the angular mesh.
        coherent : bool
            Flag indicating coherent detection mode.
        NA : float
            Numerical aperture exposed as a read-only property.
        sampling : int
            Sampling count exposed as a read-only property.
        phi_offset : float
            Azimuthal offset exposed as a read-only property.
        gamma_offset : float
            Polar offset exposed as a read-only property.
        polarization_filter : float
            Polarization filter setting exposed as a read-only property.
        rotation : float
            Rotation angle exposed as a read-only property.
        mesh : numpy.ndarray
            Underlying Fibonacci mesh of angles used by the detector.
        max_angle : float
            Maximum polar angle in the mesh (radians).
        min_angle : float
            Minimum polar angle in the mesh (radians); zero if no cache.

    )pbdoc")
        .def(py::init<std::string, size_t, double, double, double, double, double, double, bool, bool, double>(),
             py::arg("mode_number"),
             py::arg("sampling"),
             py::arg("NA"),
             py::arg("cache_NA"),
             py::arg("phi_offset"),
             py::arg("gamma_offset"),
             py::arg("polarization_filter"),
             py::arg("rotation"),
             py::arg("coherent"),
             py::arg("mean_coupling"),
             py::arg("medium_refractive_index") = 1.0,
             R"pbdoc(
                 Initialize a BindedDetector.

                 See class docstring for parameter descriptions.
             )pbdoc")

        .def_readonly(
            "mode_field",
            &Detector::mode_field,
            R"pbdoc(
                ModeField instance representing the detector's mode.
                Contains methods for field computation and projections.
            )pbdoc")
        .def_readonly(
            "mode_id",
            &Detector::mode_id,
            R"pbdoc(
                ModeID instance representing the detector's mode identifier.
                Contains family and number information for the mode.
            )pbdoc")
        .def("get_structured_scalarfield",
            &Detector::get_structured_scalarfield,
            py::arg("sampling"),
            py::return_value_policy::move,
            R"pbdoc(
                Get structured mode field data.

                Returns a NumPy array of shape (sampling, 2) containing
                (gamma, phi) pairs for the mode field.
            )pbdoc")

        .def("_cpp_get_coupling", &Detector::get_coupling,
             py::arg("scatterer"),
             R"pbdoc(
                 Compute coupling between this detector and a scatterer.

                 Parameters
                 ----------
                 scatterer : BaseScatterer
                     An instance of a PyMieSim scatterer (e.g., SPHERE or CYLINDER).

                 Returns
                 -------
                 float
                     The coupling coefficient (overlap integral) value.
             )pbdoc")

        .def_readwrite("_cpp_scalar_field", &Detector::scalar_field,
             R"pbdoc(
                 Scalar field intensity sampled on the detector mesh.

                 This NumPy array holds the detected intensity values at each
                 angular sample point.
             )pbdoc")

        .def_readwrite("coherent", &Detector::coherent,
             R"pbdoc(
                 Coherent detection mode flag.

                 True if the detector is in coherent mode; false for incoherent.
             )pbdoc")

        .def_readonly("_cpp_NA", &Detector::numerical_aperture,
             R"pbdoc(
                 Numerical aperture of the detector.

                 Controls the angular extent of collected light.
             )pbdoc")

        .def_readonly("_cpp_sampling", &Detector::sampling,
             R"pbdoc(
                 Number of samples per full angular sweep.

                 Determines resolution of the detected field.
             )pbdoc")

        .def_readonly("_cpp_phi_offset", &Detector::phi_offset,
             R"pbdoc(
                 Azimuthal offset angle for detector orientation.
             )pbdoc")

        .def_readonly("_cpp_gamma_offset", &Detector::gamma_offset,
             R"pbdoc(
                 Polar offset angle for detector orientation.
             )pbdoc")

        .def_readonly("_cpp_polarization_filter", &Detector::polarization_filter,
             R"pbdoc(
                 Polarization filter setting.

                 Specifies any linear or circular filter applied to the signal.
             )pbdoc")

        .def_readonly("_cpp_rotation", &Detector::rotation,
             R"pbdoc(
                 Rotation angle of the detector's field-of-view.
             )pbdoc")

        .def_readonly("_cpp_mesh", &Detector::fibonacci_mesh,
             R"pbdoc(
                 The Fibonacci angular mesh used for sampling.

                 Provided as a NumPy array of shape (sampling, 2) containing
                 (gamma, phi) pairs.
             )pbdoc")

        .def_readonly("_cpp_max_angle", &Detector::max_angle,
             R"pbdoc(
                 Maximum polar angle in the mesh (radians).
             )pbdoc")

        .def_readonly("_cpp_min_angle", &Detector::min_angle,
             R"pbdoc(
                 Minimum polar angle in the mesh (radians);
                 zero if `cache_NA` disabled.
             )pbdoc")
    ;
}
