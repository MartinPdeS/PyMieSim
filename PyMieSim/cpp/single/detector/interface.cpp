#include <pybind11/pybind11.h>
#include "./detector.h"
#include <fibonacci/interface.cpp>
#include <coordinates/interface.cpp>
#include <mode_field/interface.cpp>
#include <utils/numpy_interface.h>
#include <pint/pint.h>

namespace py = pybind11;

void register_detector(py::module_& module) {
    py::object ureg = get_shared_ureg();

    module.doc() = R"pbdoc(
        Photodiode binding for PyMieSim.

        Provides a `DETECTOR` class to compute coupling to scatterers
        and expose the detected field distribution.
    )pbdoc";

    register_coordinates(module);
    register_fibonacci(module);
    register_mode_field(module);

    py::class_<BaseDetector, std::shared_ptr<BaseDetector>>(module, "BASE_DETECTOR")
        .def_property_readonly(
            "numerical_aperture",
            [ureg](BaseDetector& self) {
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
            [ureg](BaseDetector& self) {
                return py::float_(self.max_angle) * ureg.attr("radian");
            },
            R"pbdoc(
            Returns the maximum angle of the detector in radians.
            This is used to determine the angular coverage of the detector.
            )pbdoc"
        )
        .def_property_readonly(
            "min_angle",
            [ureg](BaseDetector& self) {
                return py::float_(self.min_angle) * ureg.attr("radian");
            },
            R"pbdoc(
            Returns the minimum angle of the detector in radians.
            This is used to determine the angular coverage of the detector.
            )pbdoc"
        )
        .def_property(
            "scalar_field",
            [ureg](BaseDetector& self) {
                return vector_as_numpy_view(self, self.scalar_field);
            },
            [ureg](BaseDetector& self, py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.scalar_field, arr);
            },
            R"pbdoc(
                Complex far field samples as numpy.complex128, shape (N,)
                Getter returns a zero copy view tied to the detector lifetime
                Setter accepts any 1D array castable to complex128
            )pbdoc"
        )
        .def(
            "get_structured_scalarfield",
            [ureg](BaseDetector& self, size_t sampling) {
                std::vector<complex128> vector = self.get_structured_scalarfield(sampling);
                return vector_move_from_numpy(std::move(vector), {sampling, sampling});
            },
            py::arg("sampling"),
            R"pbdoc(
                NumPy array of shape (sampling, 2) with columns [gamma, phi].
            )pbdoc"
        )
        .def(
            "get_poynting_field",
            [ureg](BaseDetector& self, const BaseScatterer& scatterer, py::object distance) {
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
            [ureg](BaseDetector& self, const BaseScatterer& scatterer) {
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
            This method computes the energy flow by calculating the Poynting vector and summing it over the surface mesh of the scatterer.
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
            [ureg](BaseDetector& self) {
                return py::float_(self.phi_offset) * ureg.attr("radian");
            },
            R"pbdoc(
            Azimuthal offset angle for detector orientation.
            )pbdoc"
        )
        .def_property_readonly(
            "gamma_offset",
            [ureg](BaseDetector& self) {
                return py::float_(self.gamma_offset) * ureg.attr("radian");
            },
            R"pbdoc(
            Polar offset angle for detector orientation.
            )pbdoc"
        )
        .def_property_readonly(
            "polarization_filter",
            [ureg](BaseDetector& self) {
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
        .def(
            "get_coupling",
            [ureg](BaseDetector& self, const BaseScatterer& scatterer) {
                double coupling = self.get_coupling(scatterer);
                return (py::float_(coupling) * ureg.attr("watt")).attr("to_compact")();
            },
            py::arg("scatterer"),
            R"pbdoc(
            Compute the light coupling between the detector and a scatterer.

            The coupling quantifies the interaction between the field captured by the detector and the scattered field produced by the scatterer.

            Parameters
            ----------
            scatterer : BaseScatterer
                The scatterer object that interacts with the incident light, producing the scattered field.

            Returns
            -------
            Quantity
                The power coupling between the detector and the scatterer, expressed in watts.
            )pbdoc"
        );

    py::class_<Photodiode, BaseDetector, std::shared_ptr<Photodiode>>(module, "Photodiode",
        R"pbdoc(
            Photodiode for Lorenz-Mie scattering simulations.
        )pbdoc"
    )
        .def(
            py::init(
                [ureg](
                    py::object numerical_aperture,
                    py::object phi_offset,
                    py::object gamma_offset,
                    py::object medium,
                    py::object polarization_filter,
                    py::object cache_numerical_aperture,
                    std::size_t sampling
                ) {
                    double numerical_aperture_value = numerical_aperture.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();
                    double phi_offset_value = phi_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double gamma_offset_value = gamma_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double cache_numerical_aperture_value = cache_numerical_aperture.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();

                    double polarization_filter_value;
                    if (polarization_filter.is(py::none())) {
                        polarization_filter_value = std::nan("");
                    } else {
                        polarization_filter_value = polarization_filter.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    }

                    std::shared_ptr<BaseMedium> medium_value;
                    if (medium.is(py::none())) {
                        medium_value = std::make_shared<ConstantMedium>(1.0);
                    }
                    else if (py::isinstance<BaseMedium>(medium)) {
                        medium_value = medium.cast<std::shared_ptr<BaseMedium>>();
                    }
                    else if (py::hasattr(medium, "to")) {
                        double refractive_index_value = medium.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();
                        medium_value = std::make_shared<ConstantMedium>(refractive_index_value);
                    }
                    else {
                        try {
                            double refractive_index_value = medium.cast<double>();
                            medium_value = std::make_shared<ConstantMedium>(refractive_index_value);
                        }
                        catch (const py::cast_error&) {
                            throw std::runtime_error(
                                "medium must be None, a BaseMedium instance, or a real refractive index."
                            );
                        }
                    }

                    if (numerical_aperture_value < 0.0) {
                        throw std::runtime_error("Numerical aperture (NA) must be non-negative.");
                    }

                    if (numerical_aperture_value > 0.342) {
                        std::printf(
                            "Warning: Numerical aperture (NA) exceeds 0.342 for coherent detectors, which is not recommended for paraxial approximation to hold.\n"
                        );
                    }

                    return std::make_shared<Photodiode>(
                        sampling,
                        numerical_aperture_value,
                        cache_numerical_aperture_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        medium_value
                    );
                }
            ),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("medium") = py::none(),
            py::arg("polarization_filter") = py::float_(std::nan("")) * ureg.attr("degree"),
            py::arg("cache_numerical_aperture") = py::float_(0.0) * ureg.attr("dimensionless"),
            py::arg("sampling") = 200,
            R"pbdoc(
            Photodiode class representing a photodiode with a coherent light coupling mechanism.

            Parameters
            ----------
            numerical_aperture : Dimensionless
                Numerical aperture of the imaging system.
            gamma_offset : Angle
                Angle offset of the detector in the direction perpendicular to polarization.
            phi_offset : Angle
                Angle offset of the detector in the direction parallel to polarization.
            sampling : int
                Sampling rate of the far field distribution. Default is 200.
            polarization_filter : Optional[Angle]
                Angle of the polarization filter in front of the detector.
            cache_numerical_aperture : Optional[Dimensionless]
                Numerical aperture of the detector cache. Default is 0.
            medium : BaseMedium or float[RefractiveIndex]
                The medium in which the detector operates. Default is vacuum or air.
            )pbdoc"
        )
        .def(
            "print_properties",
            &Photodiode::print_properties,
            R"pbdoc(
            Print the properties of the Photodiode detector.
            )pbdoc"
        );

    py::class_<CoherentMode, BaseDetector, std::shared_ptr<CoherentMode>>(module, "CoherentMode",
        R"pbdoc(
            CoherentMode detector for Lorenz-Mie scattering simulations.
        )pbdoc"
    )
        .def(
            py::init(
                [ureg](
                    std::string mode_number,
                    py::object numerical_aperture,
                    py::object phi_offset,
                    py::object gamma_offset,
                    py::object rotation,
                    py::object medium,
                    py::object cache_numerical_aperture,
                    py::object polarization_filter,
                    py::object mean_coupling,
                    std::size_t sampling
                ) {
                    double numerical_aperture_value = numerical_aperture.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();
                    double phi_offset_value = phi_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double gamma_offset_value = gamma_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double rotation_value = rotation.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    double cache_numerical_aperture_value = cache_numerical_aperture.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();
                    bool mean_coupling_value = mean_coupling.cast<bool>();

                    double polarization_filter_value;
                    if (polarization_filter.is(py::none())) {
                        polarization_filter_value = std::nan("");
                    } else {
                        polarization_filter_value = polarization_filter.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    }

                    std::shared_ptr<BaseMedium> medium_value;
                    if (medium.is(py::none())) {
                        medium_value = std::make_shared<ConstantMedium>(1.0);
                    }
                    else if (py::isinstance<BaseMedium>(medium)) {
                        medium_value = medium.cast<std::shared_ptr<BaseMedium>>();
                    }
                    else if (py::hasattr(medium, "to")) {
                        double refractive_index_value = medium.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<double>();
                        medium_value = std::make_shared<ConstantMedium>(refractive_index_value);
                    }
                    else {
                        try {
                            double refractive_index_value = medium.cast<double>();
                            medium_value = std::make_shared<ConstantMedium>(refractive_index_value);
                        }
                        catch (const py::cast_error&) {
                            throw std::runtime_error(
                                "medium must be None, a BaseMedium instance, or a real refractive index."
                            );
                        }
                    }

                    if (numerical_aperture_value < 0.0) {
                        throw std::runtime_error("Numerical aperture (NA) must be non-negative.");
                    }

                    if (numerical_aperture_value > 0.342) {
                        std::printf(
                            "Warning: Numerical aperture (NA) exceeds 0.342 for coherent detectors, which is not recommended for paraxial approximation to hold.\n"
                        );
                    }

                    return std::make_shared<CoherentMode>(
                        mode_number,
                        sampling,
                        numerical_aperture_value,
                        cache_numerical_aperture_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        rotation_value,
                        mean_coupling_value,
                        medium_value
                    );
                }
            ),
            py::arg("mode_number"),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation"),
            py::arg("medium") = py::none(),
            py::arg("cache_numerical_aperture") = py::float_(0.0) * ureg.attr("dimensionless"),
            py::arg("polarization_filter") = py::float_(std::nan("")) * ureg.attr("degree"),
            py::arg("mean_coupling") = false,
            py::arg("sampling") = 200,
            R"pbdoc(
            CoherentMode class representing a photodiode with a coherent light coupling mechanism.

            Parameters
            ----------
            numerical_aperture : Dimensionless
                Numerical aperture of the imaging system.
            gamma_offset : Angle
                Angle offset of the detector in the direction perpendicular to polarization.
            phi_offset : Angle
                Angle offset of the detector in the direction parallel to polarization.
            rotation : Angle
                Rotation angle of the detector field of view.
            sampling : int
                Sampling rate of the far field distribution. Default is 200.
            polarization_filter : Optional[Angle]
                Angle of the polarization filter in front of the detector.
            cache_numerical_aperture : Optional[Dimensionless]
                Numerical aperture of the detector cache. Default is 0.
            mean_coupling : bool
                Indicates if the coupling mechanism is point wise or mean wise.
            medium : BaseMedium or float[RefractiveIndex]
                The medium in which the detector operates. Default is vacuum or air.
            )pbdoc"
        )
        .def_property_readonly(
            "rotation",
            [ureg](CoherentMode& self) {
                return py::float_(self.rotation) * ureg.attr("radian");
            },
            R"pbdoc(
                Rotation angle of the detector's field of view.
            )pbdoc"
        )
        .def(
            "print_properties",
            &CoherentMode::print_properties,
            R"pbdoc(
            Print the properties of the CoherentMode detector.
            )pbdoc"
        );

    py::class_<IntegratingSphere, BaseDetector, std::shared_ptr<IntegratingSphere>>(module, "IntegratingSphere",
        R"pbdoc(
        IntegratingSphere detector.

        4π incoherent power collector.
        )pbdoc"
    )
        .def(
            py::init(
                [ureg](
                    py::object polarization_filter,
                    std::size_t sampling
                ) {
                    double polarization_filter_value;
                    if (polarization_filter.is(py::none())) {
                        polarization_filter_value = std::nan("");
                    } else {
                        polarization_filter_value = polarization_filter.attr("to")(ureg.attr("radian")).attr("magnitude").cast<double>();
                    }

                    return std::make_shared<IntegratingSphere>(
                        sampling,
                        polarization_filter_value
                    );
                }
            ),
            py::arg("polarization_filter") = py::none(),
            py::arg("sampling") = 400,
            R"pbdoc(
            Construct an IntegratingSphere detector.

            Parameters
            ----------
            polarization_filter : Optional[Angle]
                Polarization filter angle. If None, no polarization filtering is applied.
            sampling : int
                Number of angular samples on the Fibonacci mesh.
            )pbdoc"
        )
        .def(
            "print_properties",
            &IntegratingSphere::print_properties,
            R"pbdoc(
            Print the properties of the IntegratingSphere detector.
            )pbdoc"
        );
}