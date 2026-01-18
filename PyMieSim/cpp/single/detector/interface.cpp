#include <pybind11/pybind11.h>
#include "./detector.h"
#include <fibonacci/interface.cpp>
#include <coordinates/interface.cpp>
#include <mode_field/interface.cpp>
#include <utils/numpy_interface.h>
#include <pint/pint.h>

namespace py = pybind11;

void register_detector(py::module_& module) {
    module.doc() = R"pbdoc(
        Photodiode binding for PyMieSim.

        Provides a `DETECTOR` class to compute coupling to scatterers
        and expose the detected field distribution.
    )pbdoc"
    ;

    register_coordinates(module);

    register_fibonacci(module);

    register_mode_field(module);

    // ------------------ Bindings for Photodiodes ------------------
    py::class_<BaseDetector, std::shared_ptr<BaseDetector>>(module, "BASE_DETECTOR")
        .def_property_readonly(
            "NA",
            [](BaseDetector& self) {
                return py::float_(self.numerical_aperture) * ureg.attr("dimensionless");
            },
            R"pbdoc(
            Numerical aperture of the detector.

            Controls the angular extent of collected light.
            )pbdoc"
        )
        .def_readonly(
            "_cpp_mesh",
            &BaseDetector::fibonacci_mesh,
            R"pbdoc(
            The Fibonacci angular mesh used for sampling.

            Provided as a NumPy array of shape (sampling, 2) containing
            (gamma, phi) pairs.
            )pbdoc"
        )
        .def_property_readonly(
            "max_angle",
            [](BaseDetector& self) {
                printf("Getting max_angle\n");
                return py::float_(self.max_angle) * ureg.attr("radian");
            },
            R"pbdoc(
            Returns the maximum angle of the detector in radians.
            This is used to determine the angular coverage of the detector.
            )pbdoc"
        )
        .def_property_readonly(
            "min_angle",
            [](BaseDetector& self) {
                return py::float_(self.min_angle) * ureg.attr("radian");
            },
            R"pbdoc(
            Returns the minimum angle of the detector in radians.
            This is used to determine the angular coverage of the detector.
            )pbdoc"
        )
        .def_property(
            "scalar_field",
            [](BaseDetector& self) {return vector_as_numpy_view(self, self.scalar_field);},
            [](BaseDetector& self,
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
            [](BaseDetector& self, size_t sampling) {
                std::vector<complex128> vector = self.get_structured_scalarfield(sampling);  // returns vector<double>
                return vector_move_from_numpy(std::move(vector), {sampling, sampling});
            },
            py::arg("sampling"),
            R"pbdoc(
                NumPy array of shape (sampling, 2) with columns [gamma, phi].
            )pbdoc"
        )
        .def("get_poynting_field",
            [](BaseDetector& self, const BaseScatterer& scatterer, py::object distance) {
                double distance_value = distance.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                std::vector<double> vector = self.get_poynting_field(scatterer, distance_value);
                py::array_t<double> array = vector_move_from_numpy(std::move(vector), {vector.size()});
                return (array * ureg.attr("watt/meter**2")).attr("to_compact")();
            },
            py::arg("scatterer"),
            py::arg("distance") = 1.0,
            R"pbdoc(
                Compute the Poynting vector norm, representing the energy flux density of the electromagnetic field.

                The Poynting vector describes the directional energy transfer per unit area for an electromagnetic wave. It is defined as:

                .. math::
                    \vec{S} = \epsilon_0 c^2 \, \vec{E} \times \vec{B}

                Where:

                - \( \vec{S} \): Poynting vector (W/m²)
                - \( \epsilon_0 \): Permittivity of free space (F/m)
                - \( c \): Speed of light in vacuum (m/s)
                - \( \vec{E} \): Electric field vector (V/m)
                - \( \vec{B} \): Magnetic field vector (T)

                The cross product of the electric and magnetic field vectors results in the Poynting vector, which represents the flow of electromagnetic energy in space.

                Parameters
                ----------
                scatterer : BaseScatterer
                    An instance of a PyMieSim scatterer (e.g., SPHERE or CYLINDER).
                distance : float, optional
                    Distance from the scatterer to the detector (default is 1).

                Returns
                -------
                numpy.ndarray
                    A 2D array of shape (sampling, sampling) containing the Poynting vector magnitudes.

                Notes
                -----
                The Poynting vector is computed over a 3D mesh of voxels that cover the entire solid angle of \( 4\pi \) steradians. This method calculates the local energy flux at each voxel and returns the norm, which represents the magnitude of energy flow at each point in space around the scatterer.

                The Poynting vector is fundamental in understanding how energy is transmitted through space in the form of electromagnetic waves.

                Example
                -------
                This method is used to assess the distribution of energy around a scatterer. The total energy flow can be obtained by integrating the Poynting vector over the surface enclosing the scatterer.
            )pbdoc"
        )
        .def(
            "get_energy_flow",
            [](BaseDetector& self, const BaseScatterer& scatterer) {
                double energy_flow = self.get_energy_flow(scatterer);
                return (py::float_(energy_flow) * ureg.attr("watt")).attr("to_compact")();
            },
            py::arg("scatterer"),
            R"pbdoc(
            Calculate the total energy flow (or radiated power) from the scatterer based on the Poynting vector.

            The energy flow is computed using the following relationship between the scattered energy and the incident intensity:

            .. math::
                W_a &= \sigma_{sca} \cdot I_{inc} \\[10pt]
                P &= \int_{A} I \, dA \\[10pt]
                I &= \frac{c n \epsilon_0}{2} \, |E|^2

            Where:

            - \( W_a \): Energy flow (W)
            - \( \sigma_{sca} \): Scattering cross section (m²)
            - \( I_{inc} \): Incident intensity (W/m²)
            - \( P \): Radiated power (W)
            - \( I \): Energy density (W/m²)
            - \( c \): Speed of light in vacuum (m/s)
            - \( n \): Refractive index of the surrounding medium
            - \( \epsilon_0 \): Permittivity of free space (F/m)
            - \( E \): Electric field (V/m)

            The total power is computed by integrating the intensity over the surface area of the scatterer.

            Parameters
            ----------
            scatterer : BaseScatterer
                The scatterer object, which contains information about the scattering properties of the particle, such as geometry and material.

            Returns
            -------
            Quantity
                The total energy flow (radiated power) from the scatterer, expressed in watts.

            Notes
            -----
            This method computes the energy flow by calculating the Poynting vector (which represents the directional energy flux) and summing it over the surface mesh of the scatterer. The final result is the total radiated power.
            )pbdoc"
        )
        .def_readonly(
            "sampling",
            &BaseDetector::sampling,
            R"pbdoc(
            Number of samples per full angular sweep.

            Determines resolution of the detected field.
            )pbdoc"
        )
        .def_property_readonly(
            "phi_offset",
            [](BaseDetector& self) {
                return py::float_(self.phi_offset) * ureg.attr("radian");
            },
            R"pbdoc(
            Azimuthal offset angle for detector orientation.
            )pbdoc"
        )
        .def_property_readonly(
            "gamma_offset",
            [](BaseDetector& self) {
                return py::float_(self.gamma_offset) * ureg.attr("radian");
            },
            R"pbdoc(
            Polar offset angle for detector orientation.
            )pbdoc"
        )
        .def_property_readonly(
            "polarization_filter",
            [](BaseDetector& self) {
                return py::float_(self.polarization_filter) * ureg.attr("radian");
            },
            R"pbdoc(
            Polarization filter setting.
            Specifies any linear or circular filter applied to the signal.
            )pbdoc"
        )
        .def_readonly(
            "mode_field",
            &BaseDetector::mode_field,
            R"pbdoc(
                ModeField instance representing the detector's mode.
                Contains methods for field computation and projections.
            )pbdoc"
        )
        .def_readonly(
            "mode_id",
            &BaseDetector::mode_id,
            R"pbdoc(
                ModeID instance representing the detector's mode identifier.
                Contains family and number information for the mode.
            )pbdoc"
        )
        .def_readonly(
            "is_coherent",
            &BaseDetector::is_coherent,
            R"pbdoc(
            CoherentMode detection mode flag.

            True if the detector is in coherent mode; false for incoherent.
            )pbdoc"
        )
        .def("get_coupling",
            [](BaseDetector& self, const BaseScatterer& scatterer) {
                double coupling = self.get_coupling(scatterer);

                return (py::float_(coupling) * ureg.attr("watt")).attr("to_compact")();
            },
            py::arg("scatterer"),
            R"pbdoc(
            Compute the light coupling between the detector and a scatterer.

            The coupling quantifies the interaction between the field captured by the detector and the scattered field produced by the scatterer. Mathematically, the coupling is calculated as:

            .. math::
                |\iint_{\Omega}  \Phi_{det} \, \Psi_{scat}^* \,  d \Omega|^2

            Where:

            - \( \Phi_{det} \): The capturing field of the detector, representing the sensitivity of the detector to the incoming scattered field.
            - \( \Psi_{scat} \): The scattered field produced by the scatterer.
            - \( \Omega \): The solid angle over which the integration is performed, typically covering the full \( 4\pi \) steradians around the scatterer.
            - \( d\Omega \): The differential solid angle element.

            This integral computes the overlap between the detector's sensitivity pattern and the scattered field, which is then squared to represent the power coupled into the detector.

            Parameters
            ----------
            scatterer : BaseScatterer
                The scatterer object that interacts with the incident light, producing the scattered field.

            Returns
            -------
            Quantity
                The power coupling between the detector and the scatterer, expressed in watts (W). This value represents the amount of scattered power that is captured by the detector.

            Notes
            -----
            - The method internally invokes the appropriate binding method based on the type of scatterer (e.g., Sphere, Cylinder) to calculate the coupling.
            - The coupling depends on both the geometry of the detector and the nature of the scattered field, making it essential for evaluating the efficiency of light collection in scattering experiments.

            Example
            -------
            A common use case is to evaluate how much of the scattered light from a nanoparticle is captured by a photodiode or integrating sphere. The result can be used to estimate the efficiency of light collection for scattering measurements.
            )pbdoc"
        )
        ;

    py::class_<Photodiode, BaseDetector, std::shared_ptr<Photodiode>>(module,
        "Photodiode",
            R"pbdoc(
                Photodiode for Lorenz-Mie scattering simulations.

                This class wraps a C++ Photodiode that samples the scattered field
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
                mean_coupling : bool
                    If true, uses mean coupling; if false, uses point coupling.
                medium_refractive_index : float, optional
                    Refractive index of the surrounding medium. Default is 1.0.

                Attributes
                ----------
                scalar_field : numpy.ndarray
                    The detected scalar field (intensity) sampled on the angular mesh.
                NA : float[Dimensionless]
                    Numerical aperture exposed as a read-only property.
                sampling : int
                    Sampling count exposed as a read-only property.
                phi_offset : float[Angle]
                    Azimuthal offset exposed as a read-only property.
                gamma_offset : float[Angle]
                    Polar offset exposed as a read-only property.
                polarization_filter : float[Angle]
                    Polarization filter setting exposed as a read-only property.
                mesh : numpy.ndarray
                    Underlying Fibonacci mesh of angles used by the detector.
                max_angle : float[Angle]
                    Maximum polar angle in the mesh (radians).
                min_angle : float[Angle]
                    Minimum polar angle in the mesh (radians); zero if no cache.
            )pbdoc"
        )
        .def(
            py::init(
                [](
                    py::object NA,
                    py::object phi_offset,
                    py::object gamma_offset,
                    py::object medium_refractive_index,
                    py::object polarization_filter,
                    py::object cache_NA,
                    py::object mean_coupling,
                    std::size_t sampling
                ) {
                    py::object units_refractive_index = py::module::import("PyMieSim").attr("units").attr("RefractiveIndex")();
                    py::object units_dimensionless = py::module::import("PyMieSim").attr("units").attr("Dimensionless")();
                    py::object units_angle = py::module::import("PyMieSim").attr("units").attr("Angle")();


                    NA = units_dimensionless.attr("check")(NA);

                    polarization_filter = units_angle.attr("check")(polarization_filter);

                    cache_NA = units_dimensionless.attr("check")(cache_NA);
                    phi_offset = units_angle.attr("check")(phi_offset);
                    gamma_offset = units_angle.attr("check")(gamma_offset);
                    medium_refractive_index = units_refractive_index.attr("check")(medium_refractive_index);

                    if (!medium_refractive_index.is(py::none())) {
                        polarization_filter.attr("check")(polarization_filter);
                    } else {
                        polarization_filter = py::float_(std::nan("")) * units_angle;
                    }

                    if (!cache_NA.is(py::none())) {
                        units_dimensionless.attr("check")(cache_NA);
                    } else {
                        cache_NA = py::float_(0.0) * units_dimensionless;
                    }

                    double NA_value = NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();
                    double cache_NA_value = cache_NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();
                    double phi_offset_value = phi_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double gamma_offset_value = gamma_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double polarization_filter_value = polarization_filter.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    bool mean_coupling_value = mean_coupling.cast<bool>();
                    double medium_refractive_index_value = medium_refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();

                    if (medium_refractive_index_value <= 0.0)
                        throw std::runtime_error("Medium refractive index must be positive and non-zero.");
                    if (NA_value < 0.0)
                        throw std::runtime_error("Numerical aperture (NA) must be non-negative.");
                    if (NA_value >= medium_refractive_index_value)
                        throw std::runtime_error("Numerical aperture (NA) must be less than the medium refractive index.");
                    if (NA_value > 0.342)
                        printf("Warning: Numerical aperture (NA) exceeds 0.342 for coherent detectors, which is not recommended for paraxial approximation to hold.\n");

                    return std::make_shared<Photodiode>(
                        "NC00",
                        sampling,
                        NA_value,
                        cache_NA_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        mean_coupling_value,
                        medium_refractive_index_value
                    );
                }
            ),
            py::arg("NA"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("medium_refractive_index") = py::float_(1.0) * ureg.attr("RIU"),
            py::arg("polarization_filter") = py::float_(std::nan("")) * ureg.attr("degree"),
            py::arg("cache_NA") = py::float_(0.0) * ureg.attr("dimensionless"),
            py::arg("mean_coupling") = false,
            py::arg("sampling") = 200,
            R"pbdoc(
            Photodiode class representing a photodiode with a coherent light coupling mechanism.
            This means it is dependent on the phase of the impinging scattered light field.
            Parameters
            ----------
            NA : Dimensionless
                Numerical aperture of the imaging system.
            gamma_offset : Angle
                Angle [Degree] offset of the detector in the direction perpendicular to polarization.
            phi_offset : Angle
                Angle [Degree] offset of the detector in the direction parallel to polarization.
            sampling : int
                Sampling rate of the far-field distribution. Default is 200.
            polarization_filter : Optional[Angle]
                Angle [Degree] of the polarization filter in front of the detector.
            cache_NA : Optional[Dimensionless]
                Numerical aperture of the detector cache. Default is 0 AU.
            mean_coupling : bool
                Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
            medium_refractive_index : RefractiveIndex
                The refractive index of the medium in which the detector operates. This is important for
                determining the acceptance cone of light. Default is 1.0 (vacuum or air).
            )pbdoc"
        )
    ;

    module.def(
        "IntegratingSphere",
        [](py::object medium_refractive_index, py::object polarization_filter, std::size_t sampling) {

            py::object ureg = get_shared_ureg();
            py::object units_refractive_index = py::module::import("PyMieSim").attr("units").attr("RefractiveIndex")();
            py::object units_angle = py::module::import("PyMieSim").attr("units").attr("Angle")();

            medium_refractive_index = units_refractive_index.attr("check")(medium_refractive_index);

            if (!polarization_filter.is(py::none()))
                polarization_filter = units_angle.attr("check")(polarization_filter);
            else
                polarization_filter = py::float_(std::nan("")) * units_angle;

            const double polarization_filter_value =
                polarization_filter.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();

            const double medium_refractive_index_value =
                medium_refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();


            return std::make_shared<Photodiode>(
                "NC00",
                sampling,
                2.0,
                0.0,
                0.0,
                0.0,
                polarization_filter_value,
                false,
                1.0
            );
        },
        py::arg("medium_refractive_index") = py::float_(1.0) * get_shared_ureg().attr("RIU"),
        py::arg("polarization_filter") = py::float_(std::nan("")) * get_shared_ureg().attr("degree"),
        py::arg("sampling") = 500
    );


    // py::class_<Photodiode, BaseDetector, std::shared_ptr<Photodiode>>(module,
    //     "IntegratingSphere",
    //         R"pbdoc(
    //             IntegratingSphere for Lorenz-Mie scattering simulations.

    //             This class wraps a C++ Photodiode that samples the scattered field
    //             according to numerical aperture, orientation offsets, polarization
    //             filtering, and optional coherent or mean coupling strategies.

    //             Parameters
    //             ----------
    //             sampling : int
    //                 Number of sample points on the detector's angular mesh.
    //             polarization_filter : float
    //                 Degree or type of polarization filtering applied.
    //             medium_refractive_index : float, optional
    //                 Refractive index of the surrounding medium. Default is 1.0.


    //             Attributes
    //             ----------
    //             scalar_field : numpy.ndarray
    //                 The detected scalar field (intensity) sampled on the angular mesh.
    //             sampling : int
    //                 Sampling count exposed as a read-only property.
    //             polarization_filter : float[Angle]
    //                 Polarization filter setting exposed as a read-only property.
    //             mesh : numpy.ndarray
    //                 Underlying Fibonacci mesh of angles used by the detector.
    //             medium_refractive_index : float[RefractiveIndex]
    //                 Refractive index of the surrounding medium.
    //         )pbdoc"
    //     )
    //     .def(
    //         py::init(
    //             [](
    //                 py::object medium_refractive_index,
    //                 py::object polarization_filter,
    //                 std::size_t sampling
    //             ) {
    //                 py::object units_refractive_index = py::module::import("PyMieSim").attr("units").attr("RefractiveIndex")();
    //                 py::object units_dimensionless = py::module::import("PyMieSim").attr("units").attr("Dimensionless")();
    //                 py::object units_angle = py::module::import("PyMieSim").attr("units").attr("Angle")();
    //                 medium_refractive_index = units_refractive_index.attr("check")(medium_refractive_index);

    //                 if (!polarization_filter.is(py::none()))
    //                     polarization_filter.attr("check")(polarization_filter);
    //                 else
    //                     polarization_filter = py::float_(std::nan("")) * units_angle;

    //                 double polarization_filter_value = polarization_filter.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
    //                 double medium_refractive_index_value = medium_refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();

    //                 return std::make_shared<Photodiode>(
    //                     "NC00",
    //                     sampling,
    //                     2.0,
    //                     0.0,
    //                     0.0,
    //                     0.0,
    //                     polarization_filter_value,
    //                     false,
    //                     medium_refractive_index_value
    //                 );
    //             }
    //         ),
    //         py::arg("medium_refractive_index") = py::float_(1.0) * ureg.attr("RIU"),
    //         py::arg("polarization_filter") = py::float_(std::nan("")) * ureg.attr("degree"),
    //         py::arg("sampling") = 500,
    //         R"pbdoc(
    //         CoherentMode class representing a photodiode with a coherent light coupling mechanism.
    //         This means it is dependent on the phase of the impinging scattered light field.
    //         Parameters
    //         ----------
    //         sampling : int
    //             Sampling rate of the far-field distribution. Default is 500.
    //         polarization_filter : Optional[Angle]
    //             Angle [Degree] of the polarization filter in front of the detector.
    //         medium_refractive_index : RefractiveIndex
    //             The refractive index of the medium in which the detector operates. This is important for
    //             determining the acceptance cone of light. Default is 1.0 (vacuum or air).
    //         )pbdoc"
    //     )
    // ;

    py::class_<CoherentMode, BaseDetector, std::shared_ptr<CoherentMode>>(module,
        "CoherentMode",
            R"pbdoc(
                Photodiode for Lorenz-Mie scattering simulations.

                This class wraps a C++ Photodiode that samples the scattered field
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
                mean_coupling : bool
                    If true, uses mean coupling; if false, uses point coupling.
                medium_refractive_index : float, optional
                    Refractive index of the surrounding medium. Default is 1.0.

                Attributes
                ----------
                scalar_field : numpy.ndarray
                    The detected scalar field (intensity) sampled on the angular mesh.
                NA : float[Dimensionless]
                    Numerical aperture exposed as a read-only property.
                sampling : int
                    Sampling count exposed as a read-only property.
                phi_offset : float[Angle]
                    Azimuthal offset exposed as a read-only property.
                gamma_offset : float[Angle]
                    Polar offset exposed as a read-only property.
                polarization_filter : float[Angle]
                    Polarization filter setting exposed as a read-only property.
                rotation : float[Angle]
                    Rotation angle exposed as a read-only property.
                mesh : numpy.ndarray
                    Underlying Fibonacci mesh of angles used by the detector.
                max_angle : float[Angle]
                    Maximum polar angle in the mesh (radians).
                min_angle : float[Angle]
                    Minimum polar angle in the mesh (radians); zero if no cache.
            )pbdoc"
        )
        .def(
            py::init(
                [](
                    std::string mode_number,
                    py::object NA,
                    py::object phi_offset,
                    py::object gamma_offset,
                    py::object rotation,
                    py::object medium_refractive_index,
                    py::object cache_NA,
                    py::object polarization_filter,
                    py::object mean_coupling,
                    std::size_t sampling
                ) {
                    py::object units_refractive_index = py::module::import("PyMieSim").attr("units").attr("RefractiveIndex")();
                    py::object units_dimensionless = py::module::import("PyMieSim").attr("units").attr("Dimensionless")();
                    py::object units_angle = py::module::import("PyMieSim").attr("units").attr("Angle")();


                    NA = units_dimensionless.attr("check")(NA);

                    polarization_filter = units_angle.attr("check")(polarization_filter);

                    cache_NA = units_dimensionless.attr("check")(cache_NA);
                    phi_offset = units_angle.attr("check")(phi_offset);
                    gamma_offset = units_angle.attr("check")(gamma_offset);
                    medium_refractive_index = units_refractive_index.attr("check")(medium_refractive_index);
                    rotation = units_angle.attr("check")(rotation);

                    if (!polarization_filter.is(py::none()))
                        units_angle.attr("check")(polarization_filter);
                    else
                        polarization_filter = py::float_(std::nan("")) * units_angle;

                    if (!cache_NA.is(py::none()))
                        units_dimensionless.attr("check")(cache_NA);
                    else
                        cache_NA = py::float_(0.0) * units_dimensionless;

                    double NA_value = NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();
                    double cache_NA_value = cache_NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();
                    double phi_offset_value = phi_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double gamma_offset_value = gamma_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double polarization_filter_value = polarization_filter.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double rotation_value = rotation.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    bool mean_coupling_value = mean_coupling.cast<bool>();
                    double medium_refractive_index_value = medium_refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();

                    if (medium_refractive_index_value <= 0.0)
                        throw std::runtime_error("Medium refractive index must be positive and non-zero.");
                    if (NA_value < 0.0)
                        throw std::runtime_error("Numerical aperture (NA) must be non-negative.");
                    if (NA_value >= medium_refractive_index_value)
                        throw std::runtime_error("Numerical aperture (NA) must be less than the medium refractive index.");
                    if (NA_value > 0.342)
                        printf("Warning: Numerical aperture (NA) exceeds 0.342 for coherent detectors, which is not recommended for paraxial approximation to hold.\n");

                    return std::make_shared<CoherentMode>(
                        mode_number,
                        sampling,
                        NA_value,
                        cache_NA_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        rotation_value,
                        mean_coupling_value,
                        medium_refractive_index_value
                    );
                }
            ),
            py::arg("mode_number"),
            py::arg("NA"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation"),
            py::arg("medium_refractive_index") = py::float_(1.0) * ureg.attr("RIU"),
            py::arg("cache_NA") = py::float_(0.0) * ureg.attr("dimensionless"),
            py::arg("polarization_filter") = py::float_(std::nan("")) * ureg.attr("degree"),
            py::arg("mean_coupling") = false,
            py::arg("sampling") = 200,
            R"pbdoc(
            CoherentMode class representing a photodiode with a coherent light coupling mechanism.
            This means it is dependent on the phase of the impinging scattered light field.

            Parameters
            ----------
            NA : Dimensionless
                Numerical aperture of the imaging system.
            gamma_offset : Angle
                Angle [Degree] offset of the detector in the direction perpendicular to polarization.
            phi_offset : Angle
                Angle [Degree] offset of the detector in the direction parallel to polarization.
            sampling : int
                Sampling rate of the far-field distribution. Default is 200.
            polarization_filter : Optional[Angle]
                Angle [Degree] of the polarization filter in front of the detector.
            cache_NA : Optional[Dimensionless]
                Numerical aperture of the detector cache. Default is 0 AU.
            mean_coupling : bool
                Indicates if the coupling mechanism is point-wise (True) or mean-wise (False). Default is False.
            medium_refractive_index : RefractiveIndex
                The refractive index of the medium in which the detector operates. This is important for
                determining the acceptance cone of light. Default is 1.0 (vacuum or air).
            )pbdoc"
        )
        .def_property_readonly(
            "rotation",
            [](CoherentMode& self) {
                return py::float_(self.rotation) * ureg.attr("radian");
            },
            R"pbdoc(
            Rotation angle of the detector's field-of-view.
            )pbdoc"
        )
        ;
}
