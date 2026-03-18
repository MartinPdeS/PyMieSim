#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include "./detector_set.h"
#include <experiment/material_set/material_set.h>
#include <utils/casting.h>
#include <experiment/utils.h>

#include <limits>

namespace py = pybind11;


PYBIND11_MODULE(detector_set, module) {
    py::object ureg = get_shared_ureg();

    py::class_<BaseDetectorSet, std::shared_ptr<BaseDetectorSet>>(
        module,
        "BaseDetectorSet",
        R"pdoc(
            Base class for detector set definitions.

            A detector set stores one or more detector configurations used by
            PyMieSim experiments. Concrete detector set classes define the
            detector geometry, angular sampling, optional polarization filtering,
            and the surrounding medium.
        )pdoc"
    );


    py::class_<PhotodiodeSet, BaseDetectorSet, std::shared_ptr<PhotodiodeSet>>(
        module,
        "PhotodiodeSet",
        R"pdoc(
            Set of photodiode detectors.

            A ``PhotodiodeSet`` stores one or more intensity based detectors
            defined by their numerical aperture, angular offsets, sampling,
            optional polarization filtering, and surrounding medium.
        )pdoc"
    )
        .def(
            py::init(
                [](
                    const py::object& numerical_aperture,
                    const py::object& phi_offset,
                    const py::object& gamma_offset,
                    const py::object& sampling,
                    const py::object& cache_numerical_aperture,
                    const py::object& polarization_filter,
                    const py::object& medium
                ) {
                    return std::make_shared<PhotodiodeSet>(
                        Casting::cast_py_to_vector<unsigned>(sampling),
                        Casting::cast_py_to_vector<double>(numerical_aperture),
                        Casting::cast_py_to_vector<double>(cache_numerical_aperture),
                        Casting::cast_py_to_vector<double>(phi_offset, "radian"),
                        Casting::cast_py_to_vector<double>(gamma_offset, "radian"),
                        Casting::Polarization::cast_py_to_polarization_set(polarization_filter),
                        Casting::Material::create_material_set_from_pyobject<MediumSet, double, BaseMedium, ConstantMedium>(medium, "medium"),
                        false
                    );
                }
            ),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("sampling") = py::int_(200),
            py::arg("cache_numerical_aperture") = py::float_(0.0),
            py::arg("polarization_filter") = PolarizationState(),
            py::arg("medium") = ConstantMedium(1.0),
            R"pdoc(
                Initialize a photodiode detector set.

                Parameters
                ----------
                numerical_aperture : float or array-like of float
                    Numerical aperture of each detector.
                phi_offset : Angle or array-like of Angle
                    Azimuthal offset of each detector.
                gamma_offset : Angle or array-like of Angle
                    Polar offset of each detector.
                sampling : int or array-like of int, optional
                    Number of angular sampling points used to evaluate each
                    detector. Default is ``200``.
                cache_numerical_aperture : float or array-like of float, optional
                    Numerical aperture used for cached field evaluation.
                    Default is ``0.0``.
                polarization_filter : Angle or array-like of Angle, optional
                    Polarization analyzer angle. If ``None``, no polarization
                    filter is applied.
                medium : BaseMedium, MediumSet, float, or array-like, optional
                    Surrounding medium. This can be provided either as medium
                    objects or as refractive index values. Default is
                    ``ConstantMedium(1.0)``.

                Notes
                -----
                This constructor stores inputs as provided. It does not
                broadcast scalar inputs to a target size. Use
                :meth:`build_sequential` when you want explicit broadcasting.
            )pdoc"
        )
        .def_readonly_static(
            "attributes",
            &PhotodiodeSet::attributes,
            R"pdoc(
                Names of the detector parameters stored by ``PhotodiodeSet``.
            )pdoc"
        )
        .def_readonly(
            "sampling",
            &PhotodiodeSet::sampling,
            R"pdoc(
                Sampling point count for each detector.
            )pdoc"
        )
        .def_property_readonly(
            "numerical_aperture",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.numerical_aperture);
            },
            R"pdoc(
                Numerical aperture of each detector.
            )pdoc"
        )
        .def_property_readonly(
            "cache_numerical_aperture",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.cache_numerical_aperture);
            },
            R"pdoc(
                Cached numerical aperture associated with each detector.
            )pdoc"
        )
        .def_property_readonly(
            "phi_offset",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.phi_offset) * ureg.attr("radian");
            },
            R"pdoc(
                Azimuthal offset of each detector.
            )pdoc"
        )
        .def_property_readonly(
            "gamma_offset",
            [ureg](const PhotodiodeSet& self) {
                return py::cast(self.gamma_offset) * ureg.attr("radian");
            },
            R"pdoc(
                Polar offset of each detector.
            )pdoc"
        )
        .def(
            "get_mapping",
            [ureg](const PhotodiodeSet& self) {
                py::dict mapping;
                mapping["detector:sampling"] = py::cast(self.sampling);
                mapping["detector:NA"] = py::cast(self.numerical_aperture);
                mapping["detector:cache_NA"] = py::cast(self.cache_numerical_aperture);
                mapping["detector:phi_offset"] = py::cast(self.phi_offset) * ureg.attr("radian");
                mapping["detector:gamma_offset"] = py::cast(self.gamma_offset) * ureg.attr("radian");
                mapping["detector:polarization_filter"] = get_polarizationset_representation(self.polarization_filter_set);
                mapping["detector:medium"] = get_materialset_representation<MediumSet>(self.medium);
                return mapping;
            },
            R"pdoc(
                Return detector parameters as a dictionary.

                The returned mapping is formatted for PyMieSim parameter tracking
                and dataframe generation. Keys follow the
                ``detector:<parameter_name>`` convention.
            )pdoc"
        )
        .def_static(
            "build_sequential",
            [](
                const size_t& target_size,
                const py::object& numerical_aperture,
                const py::object& phi_offset,
                const py::object& gamma_offset,
                const py::object& sampling,
                const py::object& cache_numerical_aperture,
                const py::object& polarization_filter,
                const py::object& medium
            ) {
                std::vector<unsigned> sampling_value =
                    Casting::cast_py_to_broadcasted_vector<unsigned>(
                        "sampling",
                        sampling,
                        target_size
                    );

                std::vector<double> numerical_aperture_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "numerical_aperture",
                        numerical_aperture,
                        target_size
                    );

                std::vector<double> cache_numerical_aperture_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "cache_numerical_aperture",
                        cache_numerical_aperture,
                        target_size
                    );

                std::vector<double> phi_offset_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "phi_offset",
                        phi_offset,
                        target_size,
                        "radian"
                    );

                std::vector<double> gamma_offset_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "gamma_offset",
                        gamma_offset,
                        target_size,
                        "radian"
                    );

                PolarizationSet polarization_filter_set =
                    Casting::Polarization::cast_py_to_polarization_set(
                        polarization_filter,
                        target_size
                );

                MediumSet medium_set =
                    Casting::Material::create_material_set_from_pyobject<MediumSet, double, BaseMedium, ConstantMedium>(
                        medium,
                        "medium",
                        target_size
                );

                return std::make_shared<PhotodiodeSet>(
                    sampling_value,
                    numerical_aperture_value,
                    cache_numerical_aperture_value,
                    phi_offset_value,
                    gamma_offset_value,
                    polarization_filter_set,
                    medium_set,
                    true
                );
            },
            py::arg("target_size"),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("sampling") = py::int_(200),
            py::arg("cache_numerical_aperture") = py::float_(0.0),
            py::arg("polarization_filter") = PolarizationState(),
            py::arg("medium") = ConstantMedium(1.0),
            R"pdoc(
                Build a sequential photodiode detector set.

                Scalar inputs or length one arrays are broadcast to
                ``target_size``. Any non scalar input must already have length
                ``target_size``.

                Parameters
                ----------
                target_size : int
                    Number of detectors in the sequential set.
                numerical_aperture : float or array-like of float
                    Numerical aperture of each detector.
                phi_offset : Angle or array-like of Angle
                    Azimuthal offset of each detector.
                gamma_offset : Angle or array-like of Angle
                    Polar offset of each detector.
                sampling : int or array-like of int, optional
                    Number of angular sampling points used to evaluate each
                    detector. Default is ``200``.
                cache_numerical_aperture : float or array-like of float, optional
                    Numerical aperture used for cached field evaluation.
                    Default is ``0.0``.
                polarization_filter : Angle or array-like of Angle, optional
                    Polarization analyzer angle. If ``None``, the filter is
                    disabled for all detectors.
                medium : BaseMedium, MediumSet, float, or array-like, optional
                    Surrounding medium. This can be provided either as medium
                    objects or as refractive index values. Default is
                    ``ConstantMedium(1.0)``.

                Returns
                -------
                PhotodiodeSet
                    Detector set with ``is_sequential = True``.
            )pdoc"
        )
        ;


    py::class_<CoherentModeSet, BaseDetectorSet, std::shared_ptr<CoherentModeSet>>(
        module,
        "CoherentModeSet",
        R"pdoc(
            Set of coherent mode detectors.

            A ``CoherentModeSet`` stores one or more detectors used for
            coherent mode coupling calculations. Each detector is defined by its
            mode label, numerical aperture, angular offsets, rotation, optional
            polarization filtering, and surrounding medium.
        )pdoc"
    )
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
                    const bool& mean_coupling
                ) {
                    return std::make_shared<CoherentModeSet>(
                        Casting::cast_py_to_vector<std::string>(mode_number),
                        Casting::cast_py_to_vector<unsigned>(sampling),
                        Casting::cast_py_to_vector<double>(numerical_aperture),
                        Casting::cast_py_to_vector<double>(cache_numerical_aperture),
                        Casting::cast_py_to_vector<double>(phi_offset, "radian"),
                        Casting::cast_py_to_vector<double>(gamma_offset, "radian"),
                        Casting::Polarization::cast_py_to_polarization_set(polarization_filter),
                        Casting::cast_py_to_vector<double>(rotation, "radian"),
                        Casting::Material::create_material_set_from_pyobject<MediumSet, double, BaseMedium, ConstantMedium>(medium, "medium"),
                        mean_coupling,
                        false
                    );
                }
            ),
            py::arg("mode_number"),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation"),
            py::arg("cache_numerical_aperture") = py::float_(0.0),
            py::arg("sampling") = py::int_(200),
            py::arg("polarization_filter") = PolarizationState(),
            py::arg("medium") = ConstantMedium(1.0),
            py::arg("mean_coupling") = false,
            R"pdoc(
                Initialize a coherent mode detector set.

                Parameters
                ----------
                mode_number : str or array-like of str
                    Mode label or labels defining the collected coherent mode.
                numerical_aperture : float or array-like of float
                    Numerical aperture of each detector.
                phi_offset : Angle or array-like of Angle
                    Azimuthal offset of each detector.
                gamma_offset : Angle or array-like of Angle
                    Polar offset of each detector.
                rotation : Angle or array-like of Angle
                    Rotation of the detector mode basis.
                cache_numerical_aperture : float or array-like of float, optional
                    Numerical aperture used for cached field evaluation.
                    Default is ``0.0``.
                sampling : int or array-like of int, optional
                    Number of angular sampling points used to evaluate each
                    detector. Default is ``200``.
                polarization_filter : Angle or array-like of Angle, optional
                    Polarization analyzer angle. If ``None``, no polarization
                    filter is applied.
                medium : BaseMedium, MediumSet, float, or array-like, optional
                    Surrounding medium. This can be provided either as medium
                    objects or as refractive index values. Default is
                    ``ConstantMedium(1.0)``.
                mean_coupling : bool, optional
                    If ``True``, use mean coupling rather than mode resolved
                    coupling. Default is ``False``.

                Notes
                -----
                This constructor stores inputs as provided. It does not
                broadcast scalar inputs to a target size. Use
                :meth:`build_sequential` when you want explicit broadcasting.
            )pdoc"
        )
        .def_readonly_static(
            "attributes",
            &CoherentModeSet::attributes,
            R"pdoc(
                Names of the detector parameters stored by ``CoherentModeSet``.
            )pdoc"
        )
        .def_readonly(
            "mean_coupling",
            &CoherentModeSet::mean_coupling,
            R"pdoc(
                Whether mean coupling is used.
            )pdoc"
        )
        .def_readonly(
            "coherent",
            &CoherentModeSet::coherent,
            R"pdoc(
                Flag indicating that this detector set performs coherent
                detection.
            )pdoc"
        )
        .def_readonly(
            "mode_numbers",
            &CoherentModeSet::mode_numbers,
            R"pdoc(
                Mode labels associated with each detector.
            )pdoc"
        )
        .def_readonly(
            "sampling",
            &CoherentModeSet::sampling,
            R"pdoc(
                Sampling point count for each detector.
            )pdoc"
        )
        .def_property_readonly(
            "numerical_aperture",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.numerical_aperture);
            },
            R"pdoc(
                Numerical aperture of each detector.
            )pdoc"
        )
        .def_property_readonly(
            "cache_numerical_aperture",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.cache_numerical_aperture);
            },
            R"pdoc(
                Cached numerical aperture associated with each detector.
            )pdoc"
        )
        .def_property_readonly(
            "phi_offset",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.phi_offset) * ureg.attr("radian");
            },
            R"pdoc(
                Azimuthal offset of each detector.
            )pdoc"
        )
        .def_property_readonly(
            "gamma_offset",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.gamma_offset) * ureg.attr("radian");
            },
            R"pdoc(
                Polar offset of each detector.
            )pdoc"
        )
        .def_property_readonly(
            "rotation",
            [ureg](const CoherentModeSet& self) {
                return py::cast(self.rotation) * ureg.attr("radian");
            },
            R"pdoc(
                Rotation of the detector mode basis.
            )pdoc"
        )
        .def_readonly(
            "polarization_filter",
            &CoherentModeSet::polarization_filter_set,
            R"pdoc(
                Polarization filter angle for each detector.

                Returns ``None`` when no polarization filter is defined.
            )pdoc"
        )
        .def_property_readonly(
            "medium",
            [ureg](const CoherentModeSet& self) {
                return self.medium;
            },
            R"pdoc(
                Surrounding medium associated with the detector set.
            )pdoc"
        )
        .def(
            "get_mapping",
            [ureg](const CoherentModeSet& self) {
                py::dict mapping;
                mapping["detector:mode_number"] = py::cast(self.mode_numbers);
                mapping["detector:sampling"] = py::cast(self.sampling);
                mapping["detector:NA"] = py::cast(self.numerical_aperture);
                mapping["detector:cache_NA"] = py::cast(self.cache_numerical_aperture);
                mapping["detector:phi_offset"] = py::cast(self.phi_offset) * ureg.attr("radian");
                mapping["detector:gamma_offset"] = py::cast(self.gamma_offset) * ureg.attr("radian");
                mapping["detector:polarization_filter"] = get_polarizationset_representation(self.polarization_filter_set);
                mapping["detector:rotation"] = py::cast(self.rotation) * ureg.attr("radian");
                mapping["detector:medium"] = get_materialset_representation<MediumSet>(self.medium);
                return mapping;
            },
            R"pdoc(
                Return detector parameters as a dictionary.

                The returned mapping is formatted for PyMieSim parameter tracking
                and dataframe generation. Keys follow the
                ``detector:<parameter_name>`` convention.
            )pdoc"
        )
        .def_static(
            "build_sequential",
            [ureg](
                const size_t& target_size,
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
                std::vector<std::string> mode_number_values =
                    Casting::cast_py_to_broadcasted_vector<std::string>(
                        "mode_number",
                        mode_number,
                        target_size
                    );

                std::vector<unsigned> sampling_value =
                    Casting::cast_py_to_broadcasted_vector<unsigned>(
                        "sampling",
                        sampling,
                        target_size
                    );

                std::vector<double> numerical_aperture_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "numerical_aperture",
                        numerical_aperture,
                        target_size
                    );

                std::vector<double> cache_numerical_aperture_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "cache_numerical_aperture",
                        cache_numerical_aperture,
                        target_size
                    );

                std::vector<double> phi_offset_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "phi_offset",
                        phi_offset,
                        target_size,
                        "radian"
                    );

                std::vector<double> gamma_offset_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "gamma_offset",
                        gamma_offset,
                        target_size,
                        "radian"
                    );

                std::vector<double> rotation_value =
                    Casting::cast_py_to_broadcasted_vector<double>(
                        "rotation",
                        rotation,
                        target_size,
                        "radian"
                    );

                PolarizationSet polarization_filter_set =
                    Casting::Polarization::cast_py_to_polarization_set(
                        polarization_filter,
                        target_size
                );

                MediumSet medium_set =
                    Casting::Material::create_material_set_from_pyobject<MediumSet, double, BaseMedium, ConstantMedium>(
                        medium,
                        "medium",
                        target_size
                );

                return std::make_shared<CoherentModeSet>(
                    mode_number_values,
                    sampling_value,
                    numerical_aperture_value,
                    cache_numerical_aperture_value,
                    phi_offset_value,
                    gamma_offset_value,
                    polarization_filter_set,
                    rotation_value,
                    medium_set,
                    mean_coupling,
                    true
                );
            },
            py::arg("target_size"),
            py::arg("mode_number"),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation"),
            py::arg("cache_numerical_aperture") = py::float_(0.0),
            py::arg("sampling") = py::int_(200),
            py::arg("polarization_filter") = PolarizationState(),
            py::arg("medium") = ConstantMedium(1.0),
            py::arg("mean_coupling") = false,
            R"pdoc(
                Build a sequential coherent mode detector set.

                Scalar inputs or length one arrays are broadcast to
                ``target_size``. Any non scalar input must already have length
                ``target_size``.

                Parameters
                ----------
                target_size : int
                    Number of detectors in the sequential set.
                mode_number : str or array-like of str
                    Mode label or labels defining the collected coherent mode.
                numerical_aperture : float or array-like of float
                    Numerical aperture of each detector.
                phi_offset : Angle or array-like of Angle
                    Azimuthal offset of each detector.
                gamma_offset : Angle or array-like of Angle
                    Polar offset of each detector.
                rotation : Angle or array-like of Angle
                    Rotation of the detector mode basis.
                cache_numerical_aperture : float or array-like of float, optional
                    Numerical aperture used for cached field evaluation.
                    Default is ``0.0``.
                sampling : int or array-like of int, optional
                    Number of angular sampling points used to evaluate each
                    detector. Default is ``200``.
                polarization_filter : Angle or array-like of Angle, optional
                    Polarization analyzer angle. If ``None``, the filter is
                    disabled for all detectors.
                medium : BaseMedium, MediumSet, float, or array-like, optional
                    Surrounding medium. This can be provided either as medium
                    objects or as refractive index values. Default is
                    ``ConstantMedium(1.0)``.
                mean_coupling : bool, optional
                    If ``True``, use mean coupling rather than mode resolved
                    coupling. Default is ``False``.

                Returns
                -------
                CoherentModeSet
                    Detector set with ``is_sequential = True``.
            )pdoc"
        )
        ;
}