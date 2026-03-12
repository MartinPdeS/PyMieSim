#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h>

#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include <experiment/source_set/source_set.h>
#include <experiment/sequential_broadcast.h>
#include "./base_set.h"
#include "./sphere_set.h"
#include "./cylinder_set.h"
#include "./core_shell_set.h"


namespace py = pybind11;




MaterialSet create_material_set_from_pyobject(const py::object& material, std::string material_name) {
    if (py::isinstance<MaterialSet>(material)) {
        return py::cast<MaterialSet>(material);
    }
    else if (py::hasattr(material, "units")) {
        return MaterialSet(py::cast<std::vector<complex128>>(material));
    }
    else if (py::isinstance<py::sequence>(material) && !py::isinstance<py::str>(material)) {
        std::vector<std::shared_ptr<BaseMaterial>> material_vector =
            py::cast<std::vector<std::shared_ptr<BaseMaterial>>>(material);

        return MaterialSet(material_vector);
    }
    else {
        throw std::runtime_error(
            "Invalid type for " + material_name + ". Expected MaterialSet, a quantity array, "
            "or a sequence of Material instances."
        );
    }
}

MediumSet create_medium_set_from_pyobject(const py::object& material, std::string material_name) {
    if (py::isinstance<MediumSet>(material)) {
        return py::cast<MediumSet>(material);
    }
    else if (py::hasattr(material, "units")) {
        return MediumSet(py::cast<std::vector<double>>(material));
    }
    else if (py::isinstance<py::sequence>(material) && !py::isinstance<py::str>(material)) {
        std::vector<std::shared_ptr<BaseMedium>> medium_vector =
            py::cast<std::vector<std::shared_ptr<BaseMedium>>>(material);

        return MediumSet(medium_vector);
    }
    else {
        throw std::runtime_error(
            "Invalid type for " + material_name + ". Expected MediumSet, a quantity array, "
            "or a sequence of Medium instances."
        );
    }
}







PYBIND11_MODULE(scatterer_set, module) {
    py::object ureg = get_shared_ureg();


    py::class_<ScattererSet, std::shared_ptr<ScattererSet>>(module, "BaseScattererSet")
    .def_readonly("is_sequential", &ScattererSet::is_sequential)
    ;

    py::class_<SphereSet, ScattererSet, std::shared_ptr<SphereSet>>(module, "SphereSet")
        .def(
            py::init(
                [ureg](
                    const py::object& diameter,
                    const py::object& material,
                    const py::object& medium
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    return std::make_shared<SphereSet>(
                        diameter_value,
                        create_material_set_from_pyobject(material, "material"),
                        create_medium_set_from_pyobject(medium, "medium"),
                        false
                    );
                }
            ),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium"),
            R"pdoc(
                A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their refractive_index to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                diameter : Length
                    Diameter(s) of the spherical scatterers in meters.
                material : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the spherical scatterers themselves.
                medium : List[BaseMaterial] | List[RefractiveIndex], optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def_readonly_static(
            "attributes",
            &SphereSet::attributes,
            R"pdoc(
                A list of attribute names corresponding to the properties of the SphereSet class.
                This includes 'diameter', 'material', and 'medium', which represent the key parameters
                for defining a spherical scatterer in PyMieSim simulations.
            )pdoc"
        )
        .def_readonly_static(
            "available_measure_list",
            &SphereSet::available_measure_list,
            R"pdoc(
                A list of available measures that can be calculated for the SphereSet scatterer.
                This includes various scattering and extinction efficiencies (e.g., 'Qsca', 'Qext', 'Qabs'),
                cross-sections (e.g., 'Csca', 'Cext', 'Cabs'), Mie coefficients (e.g., 'a1', 'b1'), and other
                relevant optical properties that can be derived from the spherical scatterer configuration.
            )pdoc"
        )
        .def(
            "get_mapping",
            [ureg](const SphereSet& self) {
                py::dict mapping;
                mapping["scatterer:diameter"] = py::cast(self.diameter) * ureg.attr("meter");
                mapping["scatterer:material"] = self.material;
                mapping["scatterer:medium"] = self.medium;
                return mapping;
            },
            R"pdoc(
                Generates a mapping of scatterer attributes to their corresponding values for the SphereSet instance.
                The mapping includes keys such as 'scatterer:diameter', 'scatterer:material', and 'scatterer:medium'.
            )pdoc"
        )
        .def_property_readonly(
            "diameter",
            [ureg](const SphereSet& self) {
                return py::cast(self.diameter) * ureg.attr("meter");
            }
        )
        .def_readonly(
            "material",
            &SphereSet::material
        )
        .def_readonly(
            "medium",
            &SphereSet::medium
        )
        // .def_static(
        //     "build_sequential",
        //     [ureg](
        //         const size_t& target_size,
        //         const py::object& diameter,
        //         const py::object& refractive_index,
        //         const py::object& medium_refractive_index
        //     ) {
        //         std::vector<double> diameter_value =
        //             cast_scalar_or_array_to_vector_double(diameter.attr("to")("meter").attr("magnitude"));

        //         std::vector<complex128> refractive_index_value =
        //             cast_scalar_or_array_to_vector_complex128(refractive_index.attr("to")("RIU").attr("magnitude"));

        //         std::vector<double> medium_refractive_index_value =
        //             cast_scalar_or_array_to_vector_double(medium_refractive_index.attr("to")("RIU").attr("magnitude"));

        //         diameter_value = broadcast_vector_double("diameter", diameter_value, target_size);
        //         refractive_index_value = broadcast_vector_complex128("refractive_index", refractive_index_value, target_size);
        //         medium_refractive_index_value = broadcast_vector_double("medium_refractive_index", medium_refractive_index_value, target_size);

        //         ScattererProperties scatterer_properties(refractive_index_value);
        //         MediumProperties medium(medium_refractive_index_value);

        //         return std::make_shared<SphereSet>(
        //             diameter_value,
        //             scatterer_properties,
        //             medium,
        //             true
        //         );
        //     },
        //     py::arg("total_size"),
        //     py::arg("diameter"),
        //     py::arg("refractive_index"),
        //     py::arg("medium_refractive_index"),
        //     R"pdoc(
        //         Construct a sequential SphereSet by broadcasting scalar or single element inputs.

        //         Parameters
        //         ----------
        //         total_size : int or None
        //             Target size for broadcasting. If None, uses the size of the diameter vector after conversion.
        //         diameter : Quantity or array-like Quantity
        //             Diameter(s) of the spheres. Scalars or length 1 arrays broadcast to total_size.
        //         refractive_index : List[BaseMaterial] | List[RefractiveIndex]
        //             Refractive index or indices of the spherical scatterers themselves.
        //         medium_refractive_index : List, optional
        //             BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

        //         Returns
        //         -------
        //         SphereSet
        //             Instance with is_sequential = True.
        //     )pdoc"
        // )
        ;
















    // Binding for CYLINDER::Set
    py::class_<InfiniteCylinderSet, ScattererSet, std::shared_ptr<InfiniteCylinderSet>>(module, "InfiniteCylinderSet")
        .def(
            py::init(
                [ureg](
                    const py::object& diameter,
                    const py::object& material,
                    const py::object& medium
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    return std::make_shared<InfiniteCylinderSet>(
                        diameter_value,
                        create_material_set_from_pyobject(material, "material"),
                        create_medium_set_from_pyobject(medium, "medium"),
                        false
                    );
                }
            ),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium"),
            R"pdoc(
                A data class that represents an infinite cylinder scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding infinite cylinder scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding infinite cylinder-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                diameter : Length
                    Diameter(s) of the infinite cylinder scatterers in meters.
                material : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the infinite cylinder scatterers themselves.
                medium : List[BaseMaterial] | List[RefractiveIndex], optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def_readonly_static(
            "attributes",
            &InfiniteCylinderSet::attributes,
            R"pdoc(
                A list of attribute names corresponding to the properties of the InfiniteCylinderSet class.
                This includes 'diameter', 'property', and 'medium_property', which represent the key parameters
                for defining an infinite cylinder scatterer in PyMieSim simulations.
            )pdoc"
        )
        .def_readonly_static(
            "available_measure_list",
            &InfiniteCylinderSet::available_measure_list,
            R"pdoc(
                A list of available measurement names corresponding to the properties of the InfiniteCylinderSet class.
                This includes 'Qsca', 'Qext', 'Qabs', 'Qratio', 'Qforward', 'Qback', 'Qpr', 'Csca', 'Cext', 'Cabs',
                'Cratio', 'Cforward', 'Cback', 'Cpr', 'a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'g', and 'coupling',
                which represent the key parameters for defining an infinite cylinder scatterer in PyMieSim simulations.
            )pdoc"
        )
        .def(
            "get_mapping",
            [ureg](const InfiniteCylinderSet& self) {
                py::dict mapping;
                mapping["scatterer:diameter"] = py::cast(self.diameter) * ureg.attr("meter");
                mapping["scatterer:material"] = self.material;
                mapping["scatterer:medium"] = self.medium;
                return mapping;
            },
            R"pdoc(
                Generates a mapping of scatterer attributes to their corresponding values for the InfiniteCylinderSet instance.
                The mapping includes keys such as 'scatterer:diameter', 'scatterer:refractive_index', and 'scatterer:medium_refractive_index'.
            )pdoc"
        )
        .def_property_readonly(
            "diameter",
            [ureg](const InfiniteCylinderSet& self) {
                return py::cast(self.diameter) * ureg.attr("meter");
            }
        )
        .def_readonly(
            "material",
            &InfiniteCylinderSet::material
        )
        .def_readonly(
            "medium",
            &InfiniteCylinderSet::medium
        )
    //     .def_static(
    //         "build_sequential",
    //         [ureg](
    //             const size_t& target_size,
    //             const py::object& diameter,
    //             const py::object& refractive_index,
    //             const py::object& medium_refractive_index
    //         ) {
    //             std::vector<double> diameter_value =
    //                 cast_scalar_or_array_to_vector_double(diameter.attr("to")("meter").attr("magnitude"));

    //             std::vector<complex128> refractive_index_value =
    //                 cast_scalar_or_array_to_vector_complex128(refractive_index.attr("to")("RIU").attr("magnitude"));

    //             std::vector<double> medium_refractive_index_value =
    //                 cast_scalar_or_array_to_vector_double(medium_refractive_index.attr("to")("RIU").attr("magnitude"));

    //             diameter_value = broadcast_vector_double("diameter", diameter_value, target_size);
    //             refractive_index_value = broadcast_vector_complex128("refractive_index", refractive_index_value, target_size);
    //             medium_refractive_index_value = broadcast_vector_double("medium_refractive_index", medium_refractive_index_value, target_size);

    //             ScattererProperties scatterer_properties(refractive_index_value);
    //             MediumProperties medium(medium_refractive_index_value);

    //             return std::make_shared<InfiniteCylinderSet>(
    //                 diameter_value,
    //                 scatterer_properties,
    //                 medium,
    //                 true
    //             );
    //         },
    //         py::arg("total_size"),
    //         py::arg("diameter"),
    //         py::arg("refractive_index"),
    //         py::arg("medium_refractive_index"),
    //         R"pdoc(
    //             Construct a sequential InfiniteCylinderSet by broadcasting scalar or single element inputs.

    //             Parameters
    //             ----------
    //             total_size : int or None
    //                 Target size for broadcasting. If None, uses the size of the diameter vector after conversion.
    //             diameter : Quantity or array-like Quantity
    //                 Diameter(s) of the cylinders. Scalars or length 1 arrays broadcast to total_size.
    //             refractive_index : List[BaseMaterial] | List[RefractiveIndex]
    //                 Refractive index or indices of the cylindrical scatterers themselves.
    //             medium_refractive_index : List, optional
    //                 BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

    //             Returns
    //             -------
    //             InfiniteCylinderSet
    //                 Instance with is_sequential = True.
    //         )pdoc"
    //     )
        ;

















    // Binding for CORESHELL::Set
    py::class_<CoreShellSet, ScattererSet, std::shared_ptr<CoreShellSet>>(module, "CoreShellSet")
        .def(
            py::init(
                [ureg](
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_material,
                    const py::object& shell_material,
                    const py::object& medium
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        create_material_set_from_pyobject(core_material, "core_material"),
                        create_material_set_from_pyobject(shell_material, "shell_material"),
                        create_medium_set_from_pyobject(medium, "medium"),
                        false
                    );
                }
            ),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_material"),
            py::arg("shell_material"),
            py::arg("medium"),
            R"pdoc(
                A data class that represents an infinite cylinder scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding infinite cylinder scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding infinite cylinder-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                core_diameter : Length
                    Diameter(s) of the core of the infinite cylinder scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the infinite cylinder scatterers in meters.
                core_material : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the core of the infinite cylinder scatterers.
                shell_material : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the shell of the infinite cylinder scatterers.
                medium : List[BaseMaterial] | List[RefractiveIndex], optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def_readonly_static(
            "attributes",
            &CoreShellSet::attributes,
            R"pdoc(
                A list of attribute names corresponding to the properties of the CoreShellSet class.
                This includes 'core_diameter', 'shell_thickness', 'core_material', 'shell_material', and 'medium',
                which represent the key parameters for defining a core-shell scatterer in PyMieSim simulations.
            )pdoc"
        )
        .def_readonly_static(
            "available_measure_list",
            &CoreShellSet::available_measure_list,
            R"pdoc(
                A list of available measurement names corresponding to the properties of the CoreShellSet class.
                This includes 'Qsca', 'Qext', 'Qabs', 'Qratio', 'Qforward', 'Qback', 'Qpr', 'Csca', 'Cext', 'Cabs',
                'Cratio', 'Cforward', 'Cback', 'Cpr', 'a1', 'a2', 'a3', 'b1', 'b2', 'b3', 'g', and 'coupling',
                which represent the key parameters for defining a core-shell scatterer in PyMieSim simulations.
            )pdoc"
        )
        .def(
            "get_mapping",
            [ureg](const CoreShellSet& self) {
                py::dict mapping;
                mapping["scatterer:core_diameter"] = py::cast(self.core_diameter) * ureg.attr("meter");
                mapping["scatterer:shell_thickness"] = py::cast(self.shell_thickness) * ureg.attr("meter");
                mapping["scatterer:core_material"] = self.core_material;
                mapping["scatterer:shell_material"] = self.shell_material;
                mapping["scatterer:medium"] = self.medium;
                return mapping;
            },
            R"pdoc(
                Generates a mapping of scatterer attributes to their corresponding values for the CoreShellSet instance.
                The mapping includes keys such as 'scatterer:core_diameter', 'scatterer:shell_thickness', 'scatterer:core_refractive_index', 'scatterer:shell_refractive_index', and 'scatterer:medium_refractive_index'.
            )pdoc"
        )
        .def_property_readonly(
            "core_diameter",
            [ureg](const CoreShellSet& self) {
                return py::cast(self.core_diameter) * ureg.attr("meter");
            }
        )
        .def_property_readonly(
            "shell_thickness",
            [ureg](const CoreShellSet& self) {
                return py::cast(self.shell_thickness) * ureg.attr("meter");
            }
        )
        .def_readonly(
            "core_material",
            &CoreShellSet::core_material
        )
        .def_readonly(
            "shell_material",
            &CoreShellSet::shell_material
        )
        .def_readonly(
            "medium",
            &CoreShellSet::medium
        )
    //     .def_static(
    //         "build_sequential",
    //         [ureg](
    //             const size_t& target_size,
    //             const py::object& core_diameter,
    //             const py::object& shell_thickness,
    //             const py::object& core_refractive_index,
    //             const py::object& shell_refractive_index,
    //             const py::object& medium_refractive_index
    //         ) {
    //             std::vector<double> core_diameter_value = \
    //                 cast_scalar_or_array_to_vector_double(core_diameter.attr("to")("meter").attr("magnitude"));

    //             std::vector<double> shell_thickness_value = \
    //                 cast_scalar_or_array_to_vector_double(shell_thickness.attr("to")("meter").attr("magnitude"));

    //             std::vector<complex128> core_refractive_index_value =
    //                 cast_scalar_or_array_to_vector_complex128(core_refractive_index.attr("to")("RIU").attr("magnitude"));

    //             std::vector<complex128> shell_refractive_index_value =
    //                 cast_scalar_or_array_to_vector_complex128(shell_refractive_index.attr("to")("RIU").attr("magnitude"));

    //             std::vector<double> medium_refractive_index_value =
    //                 cast_scalar_or_array_to_vector_double(medium_refractive_index.attr("to")("RIU").attr("magnitude"));

    //             core_diameter_value = broadcast_vector_double("core_diameter", core_diameter_value, target_size);
    //             shell_thickness_value = broadcast_vector_double("shell_thickness", shell_thickness_value, target_size);
    //             core_refractive_index_value = broadcast_vector_complex128("core_refractive_index", core_refractive_index_value, target_size);
    //             shell_refractive_index_value = broadcast_vector_complex128("shell_refractive_index", shell_refractive_index_value, target_size);
    //             medium_refractive_index_value = broadcast_vector_double("medium_refractive_index", medium_refractive_index_value, target_size);

    //             return std::make_shared<CoreShellSet>(
    //                 core_diameter_value,
    //                 shell_thickness_value,
    //                 core_refractive_index_value,
    //                 shell_refractive_index_value,
    //                 medium_refractive_index_value,
    //                 true
    //             );
    //         },
    //         py::arg("total_size"),
    //         py::arg("core_diameter"),
    //         py::arg("shell_thickness"),
    //         py::arg("core_refractive_index"),
    //         py::arg("shell_refractive_index"),
    //         py::arg("medium_refractive_index"),
    //         R"pdoc(
    //             Construct a sequential CoreShellSet by broadcasting scalar or single element inputs.

    //             Parameters
    //             ----------
    //             total_size : int or None
    //                 Target size for broadcasting. If None, uses the size of the core_diameter vector after conversion.
    //             core_diameter : Quantity or array-like Quantity
    //                 Diameter(s) of the core. Scalars or length 1 arrays broadcast to total_size.
    //             shell_thickness : Quantity or array-like Quantity
    //                 Thickness(es) of the shell. Scalars or length 1 arrays broadcast to total_size.
    //             core_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
    //                 Refractive index or indices of the core.
    //             shell_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
    //                 Refractive index or indices of the shell.
    //             medium_refractive_index : List, optional
    //                 BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

    //             Returns
    //             -------
    //             CoreShellSet
    //                 Instance with is_sequential = True.
    //         )pdoc"
    //     )
        ;
}

// -
