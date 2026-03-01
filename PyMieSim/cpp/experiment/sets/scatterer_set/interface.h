#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers

#include <pint/pint.h>
#include "./scatterer_set.h"


namespace py = pybind11;


void register_scatterer_set(py::module& module) {
    py::object ureg = get_shared_ureg();


    py::class_<ScattererSet, std::shared_ptr<ScattererSet>>(module, "BaseScattererSet");

    py::class_<SphereSet, ScattererSet, std::shared_ptr<SphereSet>>(module, "Sphere")
        .def(
            py::init(
                [ureg](
                    const py::object& diameter,
                    const ScattererProperties& refractive_index,
                    const MediumProperties& medium_refractive_index,
                    bool is_sequential
                ) {

                    std::vector<double> diameter_value = \
                        cast_scalar_or_array_to_vector_double(diameter.attr("to")("meter").attr("magnitude"));

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
            py::arg("is_sequential") = false,
            "Initializes a set of spheres with given diameters, refractive indices, and medium refractive index."
        )
        .def_readonly("attributes", &SphereSet::attributes)
        .def(
            "get_mapping",
            [ureg](const SphereSet& self) {
                py::dict mapping;
                mapping["scatterer:diameter"] = py::cast(self.diameter) * ureg.attr("meter");
                mapping["scatterer:refractive_index"] = self.property;
                mapping["scatterer:medium_refractive_index"] = self.medium_property;
                return mapping;
            },
            R"pdoc(
                Generates a mapping of scatterer attributes to their corresponding values for the SphereSet instance.
                The mapping includes keys such as 'scatterer:diameter', 'scatterer:refractive_index', and 'scatterer:medium_refractive_index'.
            )pdoc"
        )
        .def_property_readonly(
            "diameter",
            [ureg](const SphereSet& self) {
                return py::cast(self.diameter) * ureg.attr("meter");
            }
        )
        .def_property_readonly(
            "refractive_index",
            [](const SphereSet& self) {
                return self.property;
            }
        )
        .def_property_readonly(
            "medium_refractive_index",
            [](const SphereSet& self) {
                return self.medium_property;
            }
        )
        ;


    // Binding for CYLINDER::Set
    py::class_<InfiniteCylinderSet, ScattererSet, std::shared_ptr<InfiniteCylinderSet>>(module, "InfiniteCylinder")
        .def(
            py::init(
                [ureg](
                    const py::object& diameter,
                    const ScattererProperties& refractive_index,
                    const MediumProperties& medium_refractive_index,
                    bool is_sequential
                ) {
                    std::vector<double> diameter_value = \
                        cast_scalar_or_array_to_vector_double(diameter.attr("to")("meter").attr("magnitude"));

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
            py::arg("is_sequential") = false,
            "Initializes a set of cylinders with given diameters, refractive indices, and medium refractive index."
        )
        .def_readonly("attributes", &InfiniteCylinderSet::attributes)
        .def(
            "get_mapping",
            [ureg](const InfiniteCylinderSet& self) {
                py::dict mapping;
                mapping["scatterer:diameter"] = py::cast(self.diameter) * ureg.attr("meter");
                mapping["scatterer:refractive_index"] = self.property;
                mapping["scatterer:medium_refractive_index"] = self.medium_property;
                return mapping;
            },
            R"pdoc(
                Generates a mapping of scatterer attributes to their corresponding values for the InfiniteCylinderSet instance.
                The mapping includes keys such as 'scatterer:diameter', 'scatterer:refractive_index', and 'scatterer:medium_refractive_index'.
            )pdoc"
        )
        ;

    // Binding for CORESHELL::Set
    py::class_<CoreShellSet, ScattererSet, std::shared_ptr<CoreShellSet>>(module, "CoreShell")
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
                    std::vector<double> core_diameter_value = \
                        cast_scalar_or_array_to_vector_double(core_diameter.attr("to")("meter").attr("magnitude"));

                    std::vector<double> shell_thickness_value = \
                        cast_scalar_or_array_to_vector_double(shell_thickness.attr("to")("meter").attr("magnitude"));


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
            py::arg("is_sequential") = false,
            "Initializes a core-shell set with specific core diameters, shell widths, core indices, shell indices, and medium refractive index."
        )
        .def_readonly("attributes", &CoreShellSet::attributes)
        .def(
            "get_mapping",
            [ureg](const CoreShellSet& self) {
                py::dict mapping;
                mapping["scatterer:core_diameter"] = py::cast(self.core_diameter);
                mapping["scatterer:shell_thickness"] = py::cast(self.shell_thickness);
                mapping["scatterer:core_refractive_index"] = self.core_property;
                mapping["scatterer:shell_refractive_index"] = self.shell_property;
                mapping["scatterer:medium_refractive_index"] = self.medium_property;
                return mapping;
            },
            R"pdoc(
                Generates a mapping of scatterer attributes to their corresponding values for the CoreShellSet instance.
                The mapping includes keys such as 'scatterer:core_diameter', 'scatterer:shell_thickness', 'scatterer:core_refractive_index', 'scatterer:shell_refractive_index', and 'scatterer:medium_refractive_index'.
            )pdoc"
        )

        ;
}

// -
