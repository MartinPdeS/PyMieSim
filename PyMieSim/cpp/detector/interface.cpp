#include <pybind11/pybind11.h>
#include "detector/detector.h"
#include "../utils/numpy_interface.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_detector, module) {
    module.doc() = R"pbdoc(
        Detector binding for PyMieSim.

        Provides a `BindedDetector` class to compute coupling to scatterers
        and expose the detected field distribution.
    )pbdoc"
    ;

    // ------------------ Bindings for VectorField ------------------
    py::class_<VectorField>(module, "VECTORFIELD")
        .def_property_readonly("data",
             [](const VectorField& self) {
                std::vector<size_t> shape = {self.sampling, 3};
                return _vector_to_numpy(self.data, shape);
            },
            R"pbdoc(
                Returns the vector field data as a NumPy array.
                This property provides access to the raw vector field data.
            )pbdoc"
        )
        .def("_cpp_get_scalar_product",
            &VectorField::get_scalar_product,
            py::arg("base"),
            R"pbdoc(
                Computes the scalar product with a base vector field.

                Parameters
                ----------
                base : VECTORFIELD
                    The base vector field to compute the scalar product with.

                Returns
                -------
                numpy.ndarray
                    A NumPy array containing the scalar products.
            )pbdoc"
        )
    ;

    // ------------------ Bindings for SphericalCoordinate ------------------
    py::class_<Cartesian>(module, "CARTESIANCOORDINATE")
        .def_property_readonly("x",
            [](const Cartesian& self) {
                return pybind11::array_t<double>(self.x.size(), self.x.data(), py::cast(self));
            },
            R"pbdoc(
                Returns x coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the x-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
        .def_property_readonly("y",
            [](const Cartesian& self) {
                return pybind11::array_t<double>(self.y.size(), self.y.data(), py::cast(self));
            },
            R"pbdoc(
                Returns y coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the y-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
        .def_property_readonly("z",
            [](const Cartesian& self) {
                return pybind11::array_t<double>(self.z.size(), self.z.data(), py::cast(self));
            },
            R"pbdoc(
                Returns z coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the z-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
    ;

    // ------------------ Bindings for SphericalCoordinate ------------------
    py::class_<Spherical>(module, "SPHERICALCOORDINATE")
        .def_property_readonly("r",
            [](const Spherical& self) {
                std::vector<size_t> shape = {self.r.size()};
                std::vector<size_t> strides = get_stride<double>(shape);

                return pybind11::array_t<double>(shape, strides, self.r.data(), py::cast(self));
            },
            R"pbdoc(
                Returns radial distances of points on the spherical mesh as a NumPy array.
                This property provides the radial distances of the mesh points in spherical coordinates.
            )pbdoc"
        )
        .def_property_readonly("phi",
            [](const Spherical& self) {
                std::vector<size_t> shape = {self.phi.size()};
                std::vector<size_t> strides = get_stride<double>(shape);

                return pybind11::array_t<double>(shape, strides, self.phi.data(), py::cast(self));
            },
            R"pbdoc(
                Returns azimuthal angles (phi) of points on the spherical mesh as a NumPy array.
                This property provides the azimuthal angles of the mesh points in spherical coordinates.
            )pbdoc"
        )
        .def_property_readonly("theta",
            [](const Spherical& self) {
                std::vector<size_t> shape = {self.theta.size()};
                std::vector<size_t> strides = get_stride<double>(shape);

                return pybind11::array_t<double>(shape, strides, self.theta.data(), py::cast(self));
            },
            R"pbdoc(
                Returns polar angles (theta) of points on the spherical mesh as a NumPy array.
                This property provides the polar angles of the mesh points in spherical coordinates.
            )pbdoc"
        )
    ;

    // ------------------ Bindings for FibonacciMesh ------------------
    py::class_<FibonacciMesh>(module, "FIBONACCIMESH")
        .def(py::init<size_t, double, double, double, double, double>(),
            py::arg("sampling"),
            py::arg("max_angle"),
            py::arg("min_angle"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation_angle"),
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
            py::arg("angle"),
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

    // ------------------ Bindings for ModeField ------------------
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
            )pbdoc"
        )
    ;

    // ------------------ Bindings for Detector ------------------
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
            )pbdoc"
        )
        .def_readonly(
            "mode_field",
            &Detector::mode_field,
            R"pbdoc(
                ModeField instance representing the detector's mode.
                Contains methods for field computation and projections.
            )pbdoc"
        )
        .def_readonly(
            "mode_id",
            &Detector::mode_id,
            R"pbdoc(
                ModeID instance representing the detector's mode identifier.
                Contains family and number information for the mode.
            )pbdoc"
        )
        .def("get_structured_scalarfield",
            &Detector::get_structured_scalarfield,
            py::arg("sampling"),
            py::return_value_policy::move,
            R"pbdoc(
                Get structured mode field data.

                Returns a NumPy array of shape (sampling, 2) containing
                (gamma, phi) pairs for the mode field.
            )pbdoc"
        )
        .def("_cpp_get_coupling",
            &Detector::get_coupling,
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
            )pbdoc"
        )
        .def_readwrite("_cpp_scalar_field",
            &Detector::scalar_field,
            R"pbdoc(
                Scalar field intensity sampled on the detector mesh.

                This NumPy array holds the detected intensity values at each
                angular sample point.
            )pbdoc"
        )
        .def_readwrite("coherent", &Detector::coherent,
            R"pbdoc(
                Coherent detection mode flag.

                True if the detector is in coherent mode; false for incoherent.
            )pbdoc"
        )
        .def_readonly("_cpp_NA", &Detector::numerical_aperture,
            R"pbdoc(
                Numerical aperture of the detector.

                Controls the angular extent of collected light.
            )pbdoc"
        )
        .def_readonly("_cpp_sampling", &Detector::sampling,
            R"pbdoc(
                Number of samples per full angular sweep.

                Determines resolution of the detected field.
            )pbdoc"
        )
        .def_readonly("_cpp_phi_offset", &Detector::phi_offset,
            R"pbdoc(
                Azimuthal offset angle for detector orientation.
            )pbdoc"
        )
        .def_readonly("_cpp_gamma_offset", &Detector::gamma_offset,
            R"pbdoc(
                Polar offset angle for detector orientation.
            )pbdoc"
        )
        .def_readonly("_cpp_polarization_filter", &Detector::polarization_filter,
            R"pbdoc(
                Polarization filter setting.

                Specifies any linear or circular filter applied to the signal.
            )pbdoc"
        )
        .def_readonly("_cpp_rotation", &Detector::rotation,
            R"pbdoc(
                Rotation angle of the detector's field-of-view.
            )pbdoc"
        )
        .def_readonly("_cpp_mesh", &Detector::fibonacci_mesh,
            R"pbdoc(
                The Fibonacci angular mesh used for sampling.

                Provided as a NumPy array of shape (sampling, 2) containing
                (gamma, phi) pairs.
            )pbdoc"
        )
        .def_readonly("_cpp_max_angle", &Detector::max_angle,
            R"pbdoc(
                Maximum polar angle in the mesh (radians).
            )pbdoc"
        )
        .def_readonly("_cpp_min_angle", &Detector::min_angle,
            R"pbdoc(
                Minimum polar angle in the mesh (radians);
                zero if `cache_NA` disabled.
            )pbdoc"
        )
    ;
}
