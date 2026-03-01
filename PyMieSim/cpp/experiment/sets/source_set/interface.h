#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <complex> // For std::complex support

#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include "./source_set.h"

using complex128 = std::complex<double>;
namespace py = pybind11;

void register_source_set(py::module& module) {
    py::object ureg = get_shared_ureg();

    py::class_<PolarizationSet, std::shared_ptr<PolarizationSet>>(module, "PolarizationSet")
        .def(
            py::init<>(
                [](const py::object& angles) {
                    std::vector<double> angles_value = \
                        cast_scalar_or_array_to_vector_double(angles.attr("to")("radian").attr("magnitude"));
                    return std::make_shared<PolarizationSet>(angles_value);
                }
            ),
            py::arg("angles")
        )
        .def(
            "__len__",
            [](const PolarizationSet& self) { return self.number_of_states(); },
            "Returns the number of polarization states defined by the angles."
        )
    ;

    // Binding for SOURCE::Set
    py::class_<BaseSourceSet, std::shared_ptr<BaseSourceSet>>(module, "BaseSourceSet");

    py::class_<GaussianSourceSet, BaseSourceSet, std::shared_ptr<GaussianSourceSet>>(module, "Gaussian")
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
                mapping["source:wavelength"] = py::cast(self.wavelength);
                mapping["source:polarization"] = py::cast(self.polarization);
                mapping["source:numerical_aperture"] = py::cast(self.numerical_aperture);
                mapping["source:optical_power"] = py::cast(self.optical_power);
                return mapping;
            },
            R"pdoc(
                Generates a mapping of source attributes to their corresponding values for the GaussianSourceSet instance.
                The mapping includes keys such as 'source:wavelength', 'source:polarization', 'source:numerical_aperture', and 'source:optical_power'.
            )pdoc"
        )
        ;

    py::class_<PlaneWaveSourceSet, BaseSourceSet, std::shared_ptr<PlaneWaveSourceSet>>(module, "PlaneWave")
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
                mapping["source:wavelength"] = py::cast(self.wavelength);
                mapping["source:polarization"] = py::cast(self.polarization);
                mapping["source:amplitude"] = py::cast(self.amplitude);
                return mapping;
            },
            R"pdoc(
                Generates a mapping of source attributes to their corresponding values for the PlaneWaveSourceSet instance.
                The mapping includes keys such as 'source:wavelength', 'source:polarization', and 'source:amplitude'.
            )pdoc"
        )
        ;
}

// -
