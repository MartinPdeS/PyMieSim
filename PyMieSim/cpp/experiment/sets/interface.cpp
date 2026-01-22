#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support

#include <pint/pint.h>
#include "properties.h"
#include "source.h"
#include "scatterer.h"
#include "detector.h"

namespace py = pybind11;
typedef std::complex<double> complex128;


void register_sets(py::module& module) {
    py::object ureg = get_shared_ureg();
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<ScattererProperties, std::shared_ptr<ScattererProperties>>(module, "ScattererProperties")
        .def(
            py::init<std::vector<complex128>>(),
            py::arg("properties")
        )
        .def(
            py::init<std::vector<std::vector<complex128>>>(),
            py::arg("properties")
        )
        ;

    py::class_<MediumProperties, std::shared_ptr<MediumProperties>>(module, "MediumProperties")
        .def(
            py::init<std::vector<double>>(),
            py::arg("properties")
        )
        .def(
            py::init<std::vector<std::vector<double>>>(),
            py::arg("properties")
        )
        ;


    py::class_<ScattererSet, std::shared_ptr<ScattererSet>>(module, "ScattererSet");

    py::class_<SphereSet, ScattererSet, std::shared_ptr<SphereSet>>(module, "SphereSet")
        .def(
            py::init(
                [ureg](
                    const py::object& diameter,
                    const ScattererProperties& refractive_index,
                    const MediumProperties& medium_refractive_index,
                    bool is_sequential
                ) {
                    py::object units_length = py::module::import("PyMieSim").attr("units").attr("Length")();

                    units_length.attr("check")(diameter);
                    std::vector<double> diameter_value = diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<SphereSet>(
                        diameter_value,
                        refractive_index,
                        medium_refractive_index,
                        is_sequential
                    );
                }
            ),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("is_sequential"),
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index."
        );


    // Binding for CYLINDER::Set
    py::class_<InfiniteCylinderSet, ScattererSet, std::shared_ptr<InfiniteCylinderSet>>(module, "InfiniteCylinderSet")
        .def(
            py::init(
                [ureg](
                    const py::object& diameter,
                    const ScattererProperties& refractive_index,
                    const MediumProperties& medium_refractive_index,
                    bool is_sequential
                ) {
                    py::object units_length = py::module::import("PyMieSim").attr("units").attr("Length")();

                    units_length.attr("check")(diameter);
                    std::vector<double> diameter_value = diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<InfiniteCylinderSet>(
                        diameter_value,
                        refractive_index,
                        medium_refractive_index,
                        is_sequential
                    );
                }
            ),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("is_sequential"),
            "Initializes a set of cylinders with given diameters, refractive indices, and medium refractive index."
        );

    // Binding for CORESHELL::Set
    py::class_<CoreShellSet, ScattererSet, std::shared_ptr<CoreShellSet>>(module, "CoreShellSet")
        .def(
            py::init(
                [ureg](
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const ScattererProperties& core_refractive_index,
                    const ScattererProperties& shell_refractive_index,
                    const MediumProperties& medium_refractive_index,
                    bool is_sequential
                ) {
                    py::object units_length = py::module::import("PyMieSim").attr("units").attr("Length")();

                    units_length.attr("check")(core_diameter);
                    std::vector<double> core_diameter_value = core_diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<std::vector<double>>();

                    units_length.attr("check")(shell_thickness);
                    std::vector<double> shell_thickness_value = shell_thickness.attr("to")(ureg.attr("meter")).attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        core_refractive_index,
                        shell_refractive_index,
                        medium_refractive_index,
                        is_sequential
                    );
                }
            ),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("is_sequential"),
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index."
        );

    // Binding for SOURCE::Set
    py::class_<BaseSourceSet, std::shared_ptr<BaseSourceSet>>(module, "BaseSourceSet");

    py::class_<GaussianSourceSet, BaseSourceSet, std::shared_ptr<GaussianSourceSet>>(module, "GaussianSourceSet")
        .def(
            py::init(
                [ureg](
                    const py::object& wavelength,
                    const std::vector<std::vector<std::complex<double>>>& jones_vector,
                    const py::object& NA,
                    const py::object& optical_power,
                    bool is_sequential
                ) {
                    py::object units_length = py::module::import("PyMieSim").attr("units").attr("Length")();
                    py::object units_dimensionless = py::module::import("PyMieSim").attr("units").attr("Dimensionless")();
                    py::object units_power = py::module::import("PyMieSim").attr("units").attr("Power")();

                    units_length.attr("check")(wavelength);
                    std::vector<double> wavelength_value = wavelength.attr("to")(ureg.attr("meter")).attr("magnitude").cast<std::vector<double>>();

                    units_dimensionless.attr("check")(NA);
                    std::vector<double> NA_values = NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_power.attr("check")(optical_power);
                    std::vector<double> optical_power_values = optical_power.attr("to")(ureg.attr("watt")).attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<GaussianSourceSet>(
                        wavelength_value,
                        jones_vector,
                        NA_values,
                        optical_power_values,
                        is_sequential
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            py::arg("is_sequential"),
            "Initializes a gaussian source set with specific wavelengths, Jones vectors, and amplitudes."
        );

    py::class_<PlaneWaveSourceSet, BaseSourceSet, std::shared_ptr<PlaneWaveSourceSet>>(module, "PlaneWaveSourceSet")
        .def(
            py::init(
                [ureg](
                    const py::object& wavelength,
                    const std::vector<std::vector<std::complex<double>>>& jones_vector,
                    const py::object& amplitude,
                    bool is_sequential
                ) {
                    py::object units_length = py::module::import("PyMieSim").attr("units").attr("Length")();
                    py::object units_dimensionless = py::module::import("PyMieSim").attr("units").attr("Dimensionless")();
                    py::object units_amplitude = py::module::import("PyMieSim").attr("units").attr("ElectricField")();

                    units_length.attr("check")(wavelength);
                    std::vector<double> wavelength_value = wavelength.attr("to")(ureg.attr("meter")).attr("magnitude").cast<std::vector<double>>();

                    units_amplitude.attr("check")(amplitude);
                    std::vector<double> amplitude_values = amplitude.attr("to")(ureg.attr("volt/meter")).attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<PlaneWaveSourceSet>(
                        wavelength_value,
                        jones_vector,
                        amplitude_values,
                        is_sequential
                    );
                }
            ),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            py::arg("is_sequential"),
            "Initializes a planewave source set with specific wavelengths, Jones vectors, and amplitudes."
        );

    py::class_<BaseDetectorSet, std::shared_ptr<BaseDetectorSet>>(module, "BaseDetectorSet");

    py::class_<PhotodiodeSet, BaseDetectorSet, std::shared_ptr<PhotodiodeSet>>(module, "PhotodiodeSet")
        .def(py::init<>())
        .def(
            py::init(
                [ureg](
                    const py::object& sampling,
                    const py::object& NA,
                    const py::object& cache_NA,
                    const py::object& phi_offset,
                    const py::object& gamma_offset,
                    const py::object& polarization_filter,
                    const py::object& medium_refractive_index,
                    bool is_sequential
                ) {
                    py::object units_length = py::module::import("PyMieSim").attr("units").attr("Length")();
                    py::object units_dimensionless = py::module::import("PyMieSim").attr("units").attr("Dimensionless")();
                    py::object units_angle = py::module::import("PyMieSim").attr("units").attr("Angle")();

                    std::vector<unsigned> sampling_value = sampling.cast<std::vector<unsigned>>();

                    units_dimensionless.attr("check")(NA);
                    std::vector<double> NA_value = NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_dimensionless.attr("check")(cache_NA);
                    std::vector<double> cache_NA_value = cache_NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_angle.attr("check")(phi_offset);
                    std::vector<double> phi_offset_value = phi_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<std::vector<double>>();

                    units_angle.attr("check")(gamma_offset);
                    std::vector<double> gamma_offset_value = gamma_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<std::vector<double>>();

                    units_dimensionless.attr("check")(polarization_filter);
                    std::vector<double> polarization_filter_value = polarization_filter.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_dimensionless.attr("check")(medium_refractive_index);
                    std::vector<double> medium_refractive_index_value = medium_refractive_index.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<PhotodiodeSet>(
                        sampling_value,
                        NA_value,
                        cache_NA_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        medium_refractive_index_value,
                        is_sequential
                    );
                }
            ),
            py::arg("sampling"),
            py::arg("NA"),
            py::arg("cache_NA"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("polarization_filter"),
            py::arg("medium_refractive_index"),
            py::arg("is_sequential"),
            "Initializes a detector set with scalar fields, numerical aperture, offsets, filters, angle, coherence, and coupling type."
        );


    py::class_<CoherentModeSet, BaseDetectorSet, std::shared_ptr<CoherentModeSet>>(module, "CoherentModeSet")
        .def(py::init<>())
        .def(
            py::init(
                [ureg](
                    const py::object& mode_number,
                    const py::object& sampling,
                    const py::object& NA,
                    const py::object& cache_NA,
                    const py::object& phi_offset,
                    const py::object& gamma_offset,
                    const py::object& polarization_filter,
                    const py::object& rotation,
                    const py::object& medium_refractive_index,
                    const bool& mean_coupling,
                    bool is_sequential
                ) {
                    py::object units_length = py::module::import("PyMieSim").attr("units").attr("Length")();
                    py::object units_dimensionless = py::module::import("PyMieSim").attr("units").attr("Dimensionless")();
                    py::object units_angle = py::module::import("PyMieSim").attr("units").attr("Angle")();

                    std::vector<std::string> mode_number_values = mode_number.cast<std::vector<std::string>>();

                    std::vector<unsigned> sampling_value = sampling.cast<std::vector<unsigned>>();

                    units_dimensionless.attr("check")(NA);
                    std::vector<double> NA_value = NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_dimensionless.attr("check")(cache_NA);
                    std::vector<double> cache_NA_value = cache_NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_angle.attr("check")(phi_offset);
                    std::vector<double> phi_offset_value = phi_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<std::vector<double>>();

                    units_angle.attr("check")(gamma_offset);
                    std::vector<double> gamma_offset_value = gamma_offset.attr("to")(ureg.attr("radian")).attr("magnitude").cast<std::vector<double>>();

                    units_dimensionless.attr("check")(polarization_filter);
                    std::vector<double> polarization_filter_value = polarization_filter.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_dimensionless.attr("check")(medium_refractive_index);
                    std::vector<double> medium_refractive_index_value = medium_refractive_index.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<std::vector<double>>();

                    units_angle.attr("check")(rotation);
                    std::vector<double> rotation_value = rotation.attr("to")(ureg.attr("radian")).attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<CoherentModeSet>(
                        mode_number_values,
                        sampling_value,
                        NA_value,
                        cache_NA_value,
                        phi_offset_value,
                        gamma_offset_value,
                        polarization_filter_value,
                        rotation_value,
                        medium_refractive_index_value,
                        mean_coupling,
                        is_sequential
                    );
                }
            ),
            py::arg("mode_number"),
            py::arg("sampling"),
            py::arg("NA"),
            py::arg("cache_NA"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("polarization_filter"),
            py::arg("rotation"),
            py::arg("medium_refractive_index"),
            py::arg("mean_coupling"),
            py::arg("is_sequential"),
            "Initializes a detector set with scalar fields, numerical aperture, offsets, filters, angle, coherence, and coupling type."
        );
}

// -
