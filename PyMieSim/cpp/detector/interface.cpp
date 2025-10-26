#include <pybind11/pybind11.h>
#include "detector/detector.h"
#include "../fibonacci/interface.cpp"
#include "../coordinates/interface.cpp"
#include "../mode_field/interface.cpp"

namespace py = pybind11;


PYBIND11_MODULE(interface_detector, module) {
    module.doc() = R"pbdoc(
        Detector binding for PyMieSim.

        Provides a `DETECTOR` class to compute coupling to scatterers
        and expose the detected field distribution.
    )pbdoc"
    ;

    register_coordinates(module);

    register_fibonacci(module);

    register_mode_field(module);


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
        is_coherent : bool
            Whether to compute coherent coupling (true) or incoherent (false).
        mean_coupling : bool
            If true, uses mean coupling; if false, uses point coupling.
        medium_refractive_index : float, optional
            Refractive index of the surrounding medium. Default is 1.0.

        Attributes
        ----------
        scalar_field : numpy.ndarray
            The detected scalar field (intensity) sampled on the angular mesh.
        is_coherent : bool
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
            py::arg("is_coherent"),
            py::arg("mean_coupling"),
            py::arg("medium_refractive_index") = 1.0,
            R"pbdoc(
                Initialize a BindedDetector.

                See class docstring for parameter descriptions.
            )pbdoc"
        )
        .def_property(
            "_cpp_scalar_field",
            [](Detector& self) {return vector_as_numpy_view(self, self.scalar_field);},
            [](Detector& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.scalar_field, arr);
            },
            R"pbdoc(
                Complex far field samples as numpy.complex128, shape (N,)
                Getter returns a zero copy view tied to the detector lifetime
                Setter accepts any 1D array castable to complex128
            )pbdoc"
        )
        .def("get_structured_scalarfield",
            [](Detector& self, size_t sampling) {
                std::vector<complex128> vector = self.get_structured_scalarfield(sampling);  // returns vector<double>
                return vector_move_from_numpy(std::move(vector), {sampling, sampling});
            },
            py::arg("sampling"),
            R"pbdoc(
                NumPy array of shape (sampling, 2) with columns [gamma, phi].
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
        .def_readwrite("is_coherent", &Detector::is_coherent,
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
