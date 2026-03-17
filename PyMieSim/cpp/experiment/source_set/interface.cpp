#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers


#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include <experiment/utils.h>
#include <experiment/source_set/source_set.h>
#include <utils/casting.h>
#include <experiment/utils.h>

namespace py = pybind11;

PYBIND11_MODULE(source_set, module) {
    py::object ureg = get_shared_ureg();

    // Binding for SOURCE::Set
    py::class_<BaseSourceSet, std::shared_ptr<BaseSourceSet>>(module, "BaseSourceSet");

    py::class_<GaussianSourceSet, BaseSourceSet, std::shared_ptr<GaussianSourceSet>>(module, "GaussianSet")
        .def(
            py::init(
                [ureg](
                    const py::object& wavelength,
                    const py::object& polarization,
                    const py::object& numerical_aperture,
                    const py::object& optical_power
                ) {
                    std::vector<double> wavelength_value = Casting::cast_py_to_vector<double>(wavelength, "meter");

                    std::vector<double> numerical_aperture_values = Casting::cast_py_to_vector<double>(numerical_aperture);

                    std::vector<double> optical_power_values = Casting::cast_py_to_vector<double>(optical_power, "watt");

                    return std::make_shared<GaussianSourceSet>(
                        wavelength_value,
                        Casting::cast_py_to_polarization_set(polarization),
                        numerical_aperture_values,
                        optical_power_values,
                        false  // is_sequential
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("numerical_aperture"),
            py::arg("optical_power"),
            R"pdoc(
                Initializes a Gaussian source set with specific wavelengths, polarization, numerical apertures, and optical powers.

                Parameters
                ----------
                wavelength : float or array-like
                    Wavelength(s) of the Gaussian source. Can be a single value or an array of values for multiple sources.
                polarization : Polarization
                    Polarization of the Gaussian source. Can be a single value or an array of values for multiple sources.
                numerical_aperture : float or array-like
                    Numerical aperture(s) of the Gaussian source. Can be a single value or an array of values for multiple sources.
                optical_power : float or array-like
                    Optical power(s) of the Gaussian source. Can be a single value or an array of values for multiple sources.
            )pdoc"
        )
        .def_readonly(
            "attributes",
            &GaussianSourceSet::attributes
        )
        .def_property_readonly(
            "wavelength",
            [ureg](const GaussianSourceSet& self) {
                return (py::cast(self.wavelength) * ureg.attr("meter"));
            },
            "Returns the wavelengths of the Gaussian source set as a NumPy array with appropriate units."
        )
        .def_property_readonly(
            "numerical_aperture",
            [ureg](const GaussianSourceSet& self) {
                return (py::cast(self.numerical_aperture));
            },
            "Returns the numerical apertures of the Gaussian source set as a NumPy array with appropriate units."
        )
        .def_property_readonly(
            "optical_power",
            [ureg](const GaussianSourceSet& self) {
                return (py::cast(self.optical_power) * ureg.attr("watt"));
            },
            "Returns the optical powers of the Gaussian source set as a NumPy array with appropriate units."
        )
        .def(
            "get_mapping",
            [ureg](const GaussianSourceSet& self) {
                py::dict mapping;
                mapping["source:wavelength"] = py::cast(self.wavelength) * ureg.attr("meter");
                mapping["source:polarization"] = get_polarizationset_representation(self.polarization);
                mapping["source:numerical_aperture"] = py::cast(self.numerical_aperture);
                mapping["source:optical_power"] = py::cast(self.optical_power) * ureg.attr("watt");
                return mapping;
            },
            R"pdoc(
                Generates a mapping of source attributes to their corresponding values for the GaussianSourceSet instance.
                The mapping includes keys such as 'source:wavelength', 'source:polarization', 'source:numerical_aperture', and 'source:optical_power'.
            )pdoc"
        )
        .def_static(
            "build_sequential",
            [ureg](
                const size_t& target_size,
                const py::object& wavelength,
                const py::object& polarization,
                const py::object& numerical_aperture,
                const py::object& optical_power
            ) {
                std::vector<double> wavelength_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "wavelength",
                        wavelength,
                        target_size,
                        "meter"
                    );

                std::vector<double> numerical_aperture_values =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "numerical_aperture",
                        numerical_aperture,
                        target_size
                    );

                std::vector<double> optical_power_values =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "optical_power",
                        optical_power,
                        target_size,
                        "watt"
                    );

                return std::make_shared<GaussianSourceSet>(
                    wavelength_value,
                    Casting::cast_py_to_polarization_set(polarization, target_size),
                    numerical_aperture_values,
                    optical_power_values,
                    true // is_sequential
                );
            },
            py::arg("target_size") = py::none(),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("numerical_aperture"),
            py::arg("optical_power"),
            R"pdoc(
                Construct a sequential GaussianSourceSet by broadcasting scalar or single element inputs.

                Parameters
                ----------
                target_size : int or None
                    Target size for broadcasting. If None, uses the maximum size among wavelength, numerical_aperture, optical_power.
                wavelength : Quantity or array-like Quantity
                    Wavelength(s). Scalars or length 1 arrays broadcast to total_size.
                polarization : PolarizationSet
                    Must have number_of_states equal to 1 or total_size.
                numerical_aperture : Quantity or array-like Quantity
                    Numerical aperture(s). Scalars or length 1 arrays broadcast to total_size.
                optical_power : Quantity or array-like Quantity
                    Optical power(s). Scalars or length 1 arrays broadcast to total_size.

                Returns
                -------
                GaussianSourceSet
                    Instance with is_sequential = True.
            )pdoc"
        )
        ;

    py::class_<PlaneWaveSourceSet, BaseSourceSet, std::shared_ptr<PlaneWaveSourceSet>>(module, "PlaneWaveSet")
        .def(
            py::init(
                [](
                    const py::object& wavelength,
                    const py::object& polarization,
                    const py::object& amplitude
                ) {
                    std::vector<double> wavelength_value = \
                        Casting::cast_py_to_vector<double>(wavelength, "meter");

                    std::vector<double> amplitude_values = \
                        Casting::cast_py_to_vector<double>(amplitude, "volt/meter");

                    return std::make_shared<PlaneWaveSourceSet>(
                        wavelength_value,
                        Casting::cast_py_to_polarization_set(polarization),
                        amplitude_values,
                        false  // is_sequential = false
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("amplitude"),
            R"pdoc(
                Initializes a plane wave source set with specific wavelengths, polarization, and amplitudes.

                Parameters
                ----------
                wavelength : float or array-like
                    Wavelength(s) of the plane wave source. Can be a single value or an array of values for multiple sources.
                polarization : Polarization
                    Polarization(s) of the plane wave source. Can be a single value or an array of values for multiple sources.
                amplitude : float or array-like
                    Amplitude(s) of the plane wave source. Can be a single value or an array of values for multiple sources.
            )pdoc"
        )
        .def_readonly(
            "attributes",
            &PlaneWaveSourceSet::attributes
        )
        .def_property_readonly(
            "wavelength",
            [ureg](const PlaneWaveSourceSet& self) {
                return (py::cast(self.wavelength) * ureg.attr("meter"));
            },
            "Returns the wavelengths of the plane wave source set as a NumPy array with appropriate units."
        )
        .def_property_readonly(
            "amplitude",
            [ureg](const PlaneWaveSourceSet& self) {
                return (py::cast(self.amplitude) * ureg.attr("volt/meter"));
            },
            "Returns the amplitudes of the plane wave source set as a NumPy array with appropriate units."
        )
        .def(
            "get_mapping",
            [ureg](const PlaneWaveSourceSet& self) {
                py::dict mapping;
                mapping["source:wavelength"] = py::cast(self.wavelength) * ureg.attr("meter");
                mapping["source:polarization"] = get_polarizationset_representation(self.polarization);
                mapping["source:amplitude"] = py::cast(self.amplitude) * ureg.attr("volt/meter");
                return mapping;
            },
            R"pdoc(
                Generates a mapping of source attributes to their corresponding values for the PlaneWaveSourceSet instance.
                The mapping includes keys such as 'source:wavelength', 'source:polarization', and 'source:amplitude'.
            )pdoc"
        )
        .def_static(
            "build_sequential",
            [ureg](
                const size_t& target_size,
                const py::object& wavelength,
                const py::object& polarization,
                const py::object& amplitude
            ) {
                std::vector<double> wavelength_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "wavelength",
                        wavelength,
                        target_size,
                        "meter"
                    );

                std::vector<double> amplitude_values =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "amplitude",
                        amplitude,
                        target_size,
                        "volt/meter"
                    );

                return std::make_shared<PlaneWaveSourceSet>(
                    wavelength_value,
                    Casting::cast_py_to_polarization_set(polarization, target_size),
                    amplitude_values,
                    true  // is_sequential = true
                );
            },
            py::arg("total_size") = py::none(),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("amplitude"),
            R"pdoc(
                Construct a sequential PlaneWaveSourceSet by broadcasting scalar or single element inputs.

                Parameters
                ----------
                total_size : int or None
                    Target size for broadcasting. If None, uses the maximum size among wavelength and amplitude.
                wavelength : Quantity or array-like Quantity
                    Wavelength(s). Scalars or length 1 arrays broadcast to total_size.
                polarization : PolarizationSet
                    Must have number_of_states equal to 1 or total_size.
                amplitude : Quantity or array-like Quantity
                    Amplitude(s). Scalars or length 1 arrays broadcast to total_size.

                Returns
                -------
                PlaneWaveSourceSet
                    Instance with is_sequential = True.
            )pdoc"
        )
        ;
}

// -
