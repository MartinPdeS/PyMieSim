#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers


#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include "./source_set.h"
#include <experiment/sequential_broadcast.h>

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
                    const std::shared_ptr<PolarizationSet>& polarization,
                    const py::object& numerical_aperture,
                    const py::object& optical_power,
                    bool is_sequential
                ) {
                    std::vector<double> wavelength_value = \
                        cast_scalar_or_array_to_vector_double(wavelength.attr("to")("meter").attr("magnitude"));

                    std::vector<double> numerical_aperture_values = \
                        cast_scalar_or_array_to_vector_double(numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                    std::vector<double> optical_power_values = \
                        cast_scalar_or_array_to_vector_double(optical_power.attr("to")("watt").attr("magnitude"));

                    return std::make_shared<GaussianSourceSet>(
                        wavelength_value,
                        polarization,
                        numerical_aperture_values,
                        optical_power_values,
                        is_sequential
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("numerical_aperture"),
            py::arg("optical_power"),
            py::arg("is_sequential") = false,
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
                is_sequential : bool, optional
                    Indicates whether the source set should be treated as sequential (default is False). If True, each index corresponds to a specific combination of parameters; if False, all combinations are generated from the provided parameter arrays.
            )pdoc"
        )
        .def_readonly("attributes", &GaussianSourceSet::attributes)
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
                return (py::cast(self.numerical_aperture) * ureg.attr("dimensionless"));
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
                mapping["source:polarization"] = py::cast(self.polarization);
                mapping["source:numerical_aperture"] = py::cast(self.numerical_aperture) * ureg.attr("dimensionless");
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
                const py::object& total_size,
                const py::object& wavelength,
                const std::shared_ptr<PolarizationSet>& polarization,
                const py::object& numerical_aperture,
                const py::object& optical_power
            ) {
                const std::optional<size_t> total_size_value = parse_optional_total_size(total_size);

                std::vector<double> wavelength_value =
                    cast_scalar_or_array_to_vector_double(wavelength.attr("to")("meter").attr("magnitude"));

                std::vector<double> numerical_aperture_values =
                    cast_scalar_or_array_to_vector_double(numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                std::vector<double> optical_power_values =
                    cast_scalar_or_array_to_vector_double(optical_power.attr("to")("watt").attr("magnitude"));

                const size_t target_size = resolve_target_size(
                    total_size_value,
                    {
                        {"wavelength", wavelength_value},
                        {"numerical_aperture", numerical_aperture_values},
                        {"optical_power", optical_power_values}
                    }
                );

                wavelength_value = broadcast_vector_double("wavelength", wavelength_value, target_size);
                numerical_aperture_values = broadcast_vector_double("numerical_aperture", numerical_aperture_values, target_size);
                optical_power_values = broadcast_vector_double("optical_power", optical_power_values, target_size);

                const size_t polarization_size = polarization ? polarization->number_of_states() : 0;
                if (polarization_size == 0) {
                    throw std::invalid_argument("polarization is null or has zero states.");
                }
                if (!(polarization_size == 1 || polarization_size == target_size)) {
                    std::ostringstream oss;
                    oss << "Inconsistent sizes: 'polarization' has number_of_states " << polarization_size
                        << " but expected 1 or " << target_size << ".";
                    throw std::invalid_argument(oss.str());
                }

                return std::make_shared<GaussianSourceSet>(
                    wavelength_value,
                    polarization,
                    numerical_aperture_values,
                    optical_power_values,
                    true
                );
            },
            py::arg("total_size") = py::none(),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("numerical_aperture"),
            py::arg("optical_power"),
            R"pdoc(
                Construct a sequential GaussianSourceSet by broadcasting scalar or single element inputs.

                Parameters
                ----------
                total_size : int or None
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
                [ureg](
                    const py::object& wavelength,
                    const std::shared_ptr<PolarizationSet>& polarization,
                    const py::object& amplitude,
                    bool is_sequential
                ) {
                    std::vector<double> wavelength_value = \
                        cast_scalar_or_array_to_vector_double(wavelength.attr("to")("meter").attr("magnitude"));

                    std::vector<double> amplitude_values = \
                        cast_scalar_or_array_to_vector_double(amplitude.attr("to")("volt/meter").attr("magnitude"));

                    return std::make_shared<PlaneWaveSourceSet>(
                        wavelength_value,
                        polarization,
                        amplitude_values,
                        is_sequential
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("amplitude"),
            py::arg("is_sequential") = false,
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
                is_sequential : bool, optional
                    Indicates whether the source set should be treated as sequential (default is False). If True, each index corresponds to a specific combination of parameters; if False, all combinations are generated from the provided parameter arrays.
            )pdoc"
        )
        .def_readonly("attributes", &PlaneWaveSourceSet::attributes)
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
                mapping["source:polarization"] = py::cast(self.polarization);
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
                const py::object& total_size,
                const py::object& wavelength,
                const std::shared_ptr<PolarizationSet>& polarization,
                const py::object& amplitude
            ) {
                const std::optional<size_t> total_size_value = parse_optional_total_size(total_size);

                std::vector<double> wavelength_value =
                    cast_scalar_or_array_to_vector_double(wavelength.attr("to")("meter").attr("magnitude"));

                std::vector<double> amplitude_values =
                    cast_scalar_or_array_to_vector_double(amplitude.attr("to")("volt/meter").attr("magnitude"));

                const size_t target_size = resolve_target_size(
                    total_size_value,
                    {
                        {"wavelength", wavelength_value},
                        {"amplitude", amplitude_values}
                    }
                );

                wavelength_value = broadcast_vector_double("wavelength", wavelength_value, target_size);
                amplitude_values = broadcast_vector_double("amplitude", amplitude_values, target_size);

                const size_t polarization_size = polarization ? polarization->number_of_states() : 0;
                if (polarization_size == 0) {
                    throw std::invalid_argument("polarization is null or has zero states.");
                }
                if (!(polarization_size == 1 || polarization_size == target_size)) {
                    std::ostringstream oss;
                    oss << "Inconsistent sizes: 'polarization' has number_of_states " << polarization_size
                        << " but expected 1 or " << target_size << ".";
                    throw std::invalid_argument(oss.str());
                }

                return std::make_shared<PlaneWaveSourceSet>(
                    wavelength_value,
                    polarization,
                    amplitude_values,
                    true
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
