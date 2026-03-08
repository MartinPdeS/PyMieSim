#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers

#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include "./detector_set.h"
#include <experiment/sequential_broadcast.h>

namespace py = pybind11;


PYBIND11_MODULE(detector_set, module) {
    py::object ureg = get_shared_ureg();

    py::class_<BaseDetectorSet, std::shared_ptr<BaseDetectorSet>>(module, "BaseDetectorSet");

    py::class_<PhotodiodeSet, BaseDetectorSet, std::shared_ptr<PhotodiodeSet>>(module, "PhotodiodeSet")
        .def(py::init<>())
        .def(
            py::init(
                [ureg](
                    const py::object& NA,
                    const py::object& phi_offset,
                    const py::object& gamma_offset,
                    const py::object& sampling,
                    const py::object& cache_NA,
                    const py::object& polarization_filter,
                    const py::object& medium,
                    bool is_sequential
                ) {
                    std::vector<unsigned> sampling_value = cast_scalar_or_array_to_vector_unsigned(sampling);

                    std::vector<double> numerical_aperture_value = cast_scalar_or_array_to_vector_double(NA.attr("to")("dimensionless").attr("magnitude"));

                    std::vector<double> cache_numerical_aperture_value = cast_scalar_or_array_to_vector_double(cache_NA.attr("to")("dimensionless").attr("magnitude"));

                    std::vector<double> phi_offset_value = cast_scalar_or_array_to_vector_double(phi_offset.attr("to")("radian").attr("magnitude"));

                    std::vector<double> gamma_offset_value = cast_scalar_or_array_to_vector_double(gamma_offset.attr("to")("radian").attr("magnitude"));

                    std::vector<double> polarization_filter_value;
                    if (polarization_filter.is_none()) {
                        polarization_filter_value = std::vector<double>(sampling_value.size(), std::nan(""));
                    } else {
                        polarization_filter_value = cast_scalar_or_array_to_vector_double(polarization_filter.attr("to")("radian").attr("magnitude"));
                    }
                    std::vector<double> medium_value = cast_scalar_or_array_to_vector_double(medium.attr("to")("RIU").attr("magnitude"));

                    return std::make_shared<PhotodiodeSet>(
                        sampling_value,
                        numerical_aperture_value,
                        cache_numerical_aperture_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        medium_value,
                        is_sequential
                    );
                }
            ),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("sampling") = py::int_(200),
            py::arg("cache_numerical_aperture") = py::float_(0.0) * ureg.attr("radian"),
            py::arg("polarization_filter") = py::none(),
            py::arg("medium") = py::float_(1.0) * ureg.attr("RIU"),
            py::arg("is_sequential") = false,
            R"pdoc(
                Initializes a detector set with scalar fields, numerical aperture, offsets, filters, angle, coherence, and coupling type.

                Parameters
                ----------
                numerical_aperture : float or array-like
                    Numerical aperture of the photodiode detector. Can be a single value or an array of values for multiple detectors.
                phi_offset : float or array-like
                    Azimuthal angle offset of the photodiode detector in radians. Can be a single value or an array of values for multiple detectors.
                gamma_offset : float or array-like
                    Polar angle offset of the photodiode detector in radians. Can be a single value or an array of values for multiple detectors.
                sampling : int or array-like, optional
                    Number of sampling points for the photodiode detector. Can be a single value or an array of values for multiple detectors. Default is 200.
                cache_numerical_aperture : float or array-like, optional
                    Numerical aperture value used for caching the photodiode detector's response. Can be a single value or an array of values for multiple detectors. Default is 0.0.          polarization_filter : float or array-like, optional
                polarization_filter : float or array-like, optional
                    Polarization filter angle for the photodiode detector in radians. Can be a single value or an array of values for multiple detectors. Default is None (no polarization filter).
                medium : float or array-like, optional
                    Refractive index of the medium in which the photodiode detector is placed. Can be a single value or an array of values for multiple detectors. Default is 1.0.
                is_sequential : bool, optional
                    Indicates whether the detector set should be treated as sequential (i.e., each detector is evaluated independently) or as a batch (i.e., all detectors are evaluated together). Default is False (batch mode).
            )pdoc"
        )
        .def_readonly_static("attributes", &PhotodiodeSet::attributes)
        .def_readonly("sampling", &PhotodiodeSet::sampling)
        .def_property_readonly(
            "numerical_aperture",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.numerical_aperture) * ureg.attr("dimensionless");
            }
        )
        .def_property_readonly(
            "cache_numerical_aperture",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.cache_numerical_aperture) * ureg.attr("dimensionless");
            }
        )
        .def_property_readonly(
            "phi_offset",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.phi_offset) * ureg.attr("radian");
            }
        )
        .def_property_readonly(
            "gamma_offset",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.gamma_offset) * ureg.attr("radian");
            }
        )
        .def(
            "get_mapping",
            [ureg](const PhotodiodeSet& self) {
                py::dict mapping;
                mapping["detector:sampling"] = py::cast(self.sampling);
                mapping["detector:NA"] = py::cast(self.numerical_aperture) * ureg.attr("dimensionless");
                mapping["detector:cache_NA"] = py::cast(self.cache_numerical_aperture) * ureg.attr("dimensionless");
                mapping["detector:phi_offset"] = py::cast(self.phi_offset) * ureg.attr("radian");
                mapping["detector:gamma_offset"] = py::cast(self.gamma_offset) * ureg.attr("radian");
                mapping["detector:polarization_filter"] = py::cast(self.polarization_filter) * ureg.attr("radian");
                mapping["detector:medium"] = py::cast(self.medium) * ureg.attr("RIU");
                return mapping;
            },
            R"pdoc(
                Generates a mapping of detector attributes to their corresponding values for the PhotodiodeSet instance.
                The mapping includes keys such as 'detector:sampling', 'detector:NA', 'detector:cache_NA', 'detector:phi_offset', 'detector:gamma_offset', 'detector:polarization_filter', and 'detector:medium'.
            )pdoc"
        )
        .def_static(
            "build_sequential",
            [ureg](
                const py::object& total_size,
                const py::object& numerical_aperture,
                const py::object& phi_offset,
                const py::object& gamma_offset,
                const py::object& sampling,
                const py::object& cache_numerical_aperture,
                const py::object& polarization_filter,
                const py::object& medium
            ) {
                const std::optional<size_t> total_size_value = parse_optional_total_size(total_size);

                std::vector<unsigned> sampling_value = cast_scalar_or_array_to_vector_unsigned(sampling);

                std::vector<double> numerical_aperture_value =
                    cast_scalar_or_array_to_vector_double(numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                std::vector<double> cache_numerical_aperture_value =
                    cast_scalar_or_array_to_vector_double(cache_numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                std::vector<double> phi_offset_value =
                    cast_scalar_or_array_to_vector_double(phi_offset.attr("to")("radian").attr("magnitude"));

                std::vector<double> gamma_offset_value =
                    cast_scalar_or_array_to_vector_double(gamma_offset.attr("to")("radian").attr("magnitude"));

                std::vector<double> medium_value =
                    cast_scalar_or_array_to_vector_double(medium.attr("to")("RIU").attr("magnitude"));

                std::vector<double> polarization_filter_value;
                if (polarization_filter.is_none()) {
                    polarization_filter_value = std::vector<double>{std::nan("")}; // broadcast later
                } else {
                    polarization_filter_value =
                        cast_scalar_or_array_to_vector_double(polarization_filter.attr("to")("radian").attr("magnitude"));
                }

                const size_t target_size = resolve_target_size(
                    total_size_value,
                    {
                        {"sampling", std::vector<double>{static_cast<double>(sampling_value.size())}}, // placeholder not used
                        {"numerical_aperture", numerical_aperture_value},
                        {"cache_numerical_aperture", cache_numerical_aperture_value},
                        {"phi_offset", phi_offset_value},
                        {"gamma_offset", gamma_offset_value},
                        {"polarization_filter", polarization_filter_value},
                        {"medium", medium_value}
                    }
                );

                sampling_value = broadcast_vector_unsigned("sampling", sampling_value, target_size);
                numerical_aperture_value = broadcast_vector_double("numerical_aperture", numerical_aperture_value, target_size);
                cache_numerical_aperture_value = broadcast_vector_double("cache_numerical_aperture", cache_numerical_aperture_value, target_size);
                phi_offset_value = broadcast_vector_double("phi_offset", phi_offset_value, target_size);
                gamma_offset_value = broadcast_vector_double("gamma_offset", gamma_offset_value, target_size);
                polarization_filter_value = broadcast_vector_double("polarization_filter", polarization_filter_value, target_size);
                medium_value = broadcast_vector_double("medium", medium_value, target_size);

                return std::make_shared<PhotodiodeSet>(
                    sampling_value,
                    numerical_aperture_value,
                    cache_numerical_aperture_value,
                    phi_offset_value,
                    gamma_offset_value,
                    polarization_filter_value,
                    medium_value,
                    true
                );
            },
            py::arg("total_size") = py::none(),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("sampling") = py::int_(200),
            py::arg("cache_numerical_aperture") = py::float_(0.0) * ureg.attr("radian"),
            py::arg("polarization_filter") = py::none(),
            py::arg("medium") = py::float_(1.0) * ureg.attr("RIU"),
            R"pdoc(
                Construct a sequential PhotodiodeSet by broadcasting scalar or single element inputs.

                total_size : int or None
                    Target size for broadcasting. If None, uses the maximum size among all provided parameters.
            )pdoc"
        )
        ;


    py::class_<CoherentModeSet, BaseDetectorSet, std::shared_ptr<CoherentModeSet>>(module, "CoherentModeSet")
        .def(py::init<>())
        .def(
            py::init(
                [ureg](
                    const py::object& mode_number,
                    const py::object& numerical_aperture,
                    const py::object& phi_offset,
                    const py::object& gamma_offset,
                    const py::object& rotation,
                    const py::object& cache_numerical_aperture,
                    const py::object& sampling,
                    const py::object& polarization_filter,
                    const py::object& medium,
                    const bool& mean_coupling,
                    bool is_sequential
                ) {
                    std::vector<std::string> mode_number_values = cast_scalar_or_array_to_vector_string(mode_number);

                    std::vector<unsigned> sampling_value = cast_scalar_or_array_to_vector_unsigned(sampling);

                    std::vector<double> numerical_aperture_value = cast_scalar_or_array_to_vector_double(numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                    std::vector<double> cache_numerical_aperture_value = cast_scalar_or_array_to_vector_double(cache_numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                    std::vector<double> phi_offset_value = cast_scalar_or_array_to_vector_double(phi_offset.attr("to")("radian").attr("magnitude"));

                    std::vector<double> gamma_offset_value = cast_scalar_or_array_to_vector_double(gamma_offset.attr("to")("radian").attr("magnitude"));

                    std::vector<double> polarization_filter_value;
                    if (polarization_filter.is_none()) {
                        polarization_filter_value = std::vector<double>(sampling_value.size(), std::nan(""));
                    } else {
                        polarization_filter_value = cast_scalar_or_array_to_vector_double(polarization_filter.attr("to")("radian").attr("magnitude"));
                    }

                    std::vector<double> medium_value = cast_scalar_or_array_to_vector_double(medium.attr("to")("RIU").attr("magnitude"));

                    std::vector<double> rotation_value = cast_scalar_or_array_to_vector_double(rotation.attr("to")("radian").attr("magnitude"));

                    return std::make_shared<CoherentModeSet>(
                        mode_number_values,
                        sampling_value,
                        numerical_aperture_value,
                        cache_numerical_aperture_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        rotation_value,
                        medium_value,
                        mean_coupling,
                        is_sequential
                    );
                }
            ),
            py::arg("mode_number"),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation"),
            py::arg("cache_numerical_aperture") = py::float_(0.0) * ureg.attr("radian"),
            py::arg("sampling") = py::int_(200),
            py::arg("polarization_filter") = py::none(),
            py::arg("medium") = py::float_(1.0) * ureg.attr("RIU"),
            py::arg("mean_coupling") = false,
            py::arg("is_sequential") = false,
            R"pdoc(
                Initializes a coherent mode detector set with mode numbers, numerical aperture, offsets, filters, angle, coherence, and coupling type.

                Parameters
                ----------
                mode_number : str or array-like
                    Identifier(s) for the coherent mode(s) to be detected. Can be a single string or an array of strings for multiple modes.
                numerical_aperture : float or array-like
                    Numerical aperture of the coherent mode detector. Can be a single value or an array of values for multiple detectors.
                phi_offset : float or array-like
                    Azimuthal angle offset of the coherent mode detector in radians. Can be a single value or an array of values for multiple detectors.
                gamma_offset : float or array-like
                    Polar angle offset of the coherent mode detector in radians. Can be a single value or an array of values for multiple detectors.
                rotation : float or array-like
                    Rotation angle of the coherent mode detector in radians. Can be a single value or an array of values for multiple detectors.
                cache_numerical_aperture : float or array-like, optional
                    Numerical aperture value used for caching the coherent mode detector's response. Can be a single value or an array of values for multiple detectors. Default is 0.0.
                sampling : int or array-like, optional
                    Number of sampling points for the coherent mode detector. Can be a single value or an array of values for multiple detectors. Default is 200.
                polarization_filter : float or array-like, optional
                    Polarization filter angle for the coherent mode detector in radians. Can be a single value or an array of values for multiple detectors. Default is NaN (no polarization filter).
                medium : float or array-like, optional
                    Refractive index of the medium in which the coherent mode detector is placed. Can be a single value or an array of values for multiple detectors. Default is 1.0.
                mean_coupling : bool, optional
                    Indicates whether to use mean coupling when calculating the response of the coherent mode detector. Default is False.
                is_sequential : bool, optional
                    Indicates whether the detector set should be treated as sequential (i.e., each detector is evaluated independently) or as a batch (i.e., all detectors are evaluated together). Default is False (batch mode).
            )pdoc"
        )
        .def_readonly_static("attributes", &CoherentModeSet::attributes)
        .def_readonly("mean_coupling", &CoherentModeSet::mean_coupling)
        .def_readonly("coherent", &CoherentModeSet::coherent)
        .def_readonly("mode_numbers", &CoherentModeSet::mode_numbers)
        .def_readonly("sampling", &CoherentModeSet::sampling)
        .def_property_readonly(
            "numerical_aperture",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.numerical_aperture) * ureg.attr("dimensionless");
            }
        )
        .def_property_readonly(
            "cache_numerical_aperture",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.cache_numerical_aperture) * ureg.attr("dimensionless");
            }
        )
        .def_property_readonly(
            "phi_offset",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.phi_offset) * ureg.attr("radian");
            }
        )
        .def_property_readonly(
            "gamma_offset",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.gamma_offset) * ureg.attr("radian");
            }
        )
        .def_property_readonly(
            "rotation",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.rotation) * ureg.attr("radian");
            }
        )
        .def_property_readonly(
            "polarization_filter",
            [ureg](const CoherentModeSet& self) {
                py::object output;
                if (std::isnan(self.polarization_filter[0])) {
                    output = py::none();
                } else {
                    output = py::cast(self.polarization_filter) * ureg.attr("dimensionless");
                }
                return output;
            }
        )
        .def_property_readonly(
            "medium",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.medium) * ureg.attr("RIU");
            }
        )
        .def(
            "get_mapping",
            [ureg](const CoherentModeSet& self) {
                py::dict mapping;
                mapping["detector:mode_number"] = py::cast(self.mode_numbers);
                mapping["detector:sampling"] = py::cast(self.sampling);
                mapping["detector:NA"] = py::cast(self.numerical_aperture) * ureg.attr("dimensionless");
                mapping["detector:cache_NA"] = py::cast(self.cache_numerical_aperture) * ureg.attr("dimensionless");
                mapping["detector:phi_offset"] = py::cast(self.phi_offset) * ureg.attr("radian");
                mapping["detector:gamma_offset"] = py::cast(self.gamma_offset) * ureg.attr("radian");
                mapping["detector:rotation"] = py::cast(self.rotation) * ureg.attr("radian");
                mapping["detector:polarization_filter"] = py::cast(self.polarization_filter) * ureg.attr("radian");
                mapping["detector:medium"] = py::cast(self.medium) * ureg.attr("RIU");
                return mapping;
            },
            R"pdoc(
                Generates a mapping of detector attributes to their corresponding values for the CoherentModeSet instance.
                The mapping includes keys such as 'detector:mode_number', 'detector:sampling', 'detector:NA', 'detector:cache_NA', 'detector:phi_offset', 'detector:gamma_offset', 'detector:polarization_filter', 'detector:rotation', 'detector:medium', and 'detector:mean_coupling'.
            )pdoc"
        )
        .def_static(
            "build_sequential",
            [ureg](
                const py::object& total_size,
                const py::object& mode_number,
                const py::object& numerical_aperture,
                const py::object& phi_offset,
                const py::object& gamma_offset,
                const py::object& rotation,
                const py::object& cache_numerical_aperture,
                const py::object& sampling,
                const py::object& polarization_filter,
                const py::object& medium,
                const bool& mean_coupling
            ) {
                const std::optional<size_t> total_size_value = parse_optional_total_size(total_size);

                std::vector<std::string> mode_number_values = cast_scalar_or_array_to_vector_string(mode_number);
                std::vector<unsigned> sampling_value = cast_scalar_or_array_to_vector_unsigned(sampling);

                std::vector<double> numerical_aperture_value =
                    cast_scalar_or_array_to_vector_double(numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                std::vector<double> cache_numerical_aperture_value =
                    cast_scalar_or_array_to_vector_double(cache_numerical_aperture.attr("to")("dimensionless").attr("magnitude"));

                std::vector<double> phi_offset_value =
                    cast_scalar_or_array_to_vector_double(phi_offset.attr("to")("radian").attr("magnitude"));

                std::vector<double> gamma_offset_value =
                    cast_scalar_or_array_to_vector_double(gamma_offset.attr("to")("radian").attr("magnitude"));

                std::vector<double> rotation_value =
                    cast_scalar_or_array_to_vector_double(rotation.attr("to")("radian").attr("magnitude"));

                std::vector<double> polarization_filter_value;
                if (polarization_filter.is_none()) {
                    polarization_filter_value = std::vector<double>{std::nan("")}; // broadcast later
                } else {
                    polarization_filter_value =
                        cast_scalar_or_array_to_vector_double(polarization_filter.attr("to")("radian").attr("magnitude"));
                }

                std::vector<double> medium_value =
                    cast_scalar_or_array_to_vector_double(medium.attr("to")("RIU").attr("magnitude"));

                const size_t target_size = resolve_target_size_from_sizes(
                    total_size_value,
                    {
                        {"mode_number", mode_number_values.size()},
                        {"sampling", sampling_value.size()},
                        {"numerical_aperture", numerical_aperture_value.size()},
                        {"cache_numerical_aperture", cache_numerical_aperture_value.size()},
                        {"phi_offset", phi_offset_value.size()},
                        {"gamma_offset", gamma_offset_value.size()},
                        {"rotation", rotation_value.size()},
                        {"polarization_filter", polarization_filter_value.size()},
                        {"medium", medium_value.size()}
                    }
                );

                mode_number_values = broadcast_vector_string("mode_number", mode_number_values, target_size);
                sampling_value = broadcast_vector_unsigned("sampling", sampling_value, target_size);
                numerical_aperture_value = broadcast_vector_double("numerical_aperture", numerical_aperture_value, target_size);
                cache_numerical_aperture_value = broadcast_vector_double("cache_numerical_aperture", cache_numerical_aperture_value, target_size);
                phi_offset_value = broadcast_vector_double("phi_offset", phi_offset_value, target_size);
                gamma_offset_value = broadcast_vector_double("gamma_offset", gamma_offset_value, target_size);
                rotation_value = broadcast_vector_double("rotation", rotation_value, target_size);
                polarization_filter_value = broadcast_vector_double("polarization_filter", polarization_filter_value, target_size);
                medium_value = broadcast_vector_double("medium", medium_value, target_size);

                return std::make_shared<CoherentModeSet>(
                    mode_number_values,
                    sampling_value,
                    numerical_aperture_value,
                    cache_numerical_aperture_value,
                    phi_offset_value,
                    gamma_offset_value,
                    polarization_filter_value,
                    rotation_value,
                    medium_value,
                    mean_coupling,
                    true
                );
            },
            py::arg("total_size") = py::none(),
            py::arg("mode_number"),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation"),
            py::arg("cache_numerical_aperture") = py::float_(0.0) * ureg.attr("radian"),
            py::arg("sampling") = py::int_(200),
            py::arg("polarization_filter") = py::none(),
            py::arg("medium") = py::float_(1.0) * ureg.attr("RIU"),
            py::arg("mean_coupling") = false,
            R"pdoc(
                Construct a sequential CoherentModeSet by broadcasting scalar or single element inputs.

                total_size : int or None
                    Target size for broadcasting. If None, uses the maximum size among all provided parameters.
            )pdoc"
        )
        ;
}

// -
