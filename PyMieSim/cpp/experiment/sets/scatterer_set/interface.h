#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers

#include <pint/pint.h>
#include "./base.h"
#include "./sphere.h"
#include "./cylinder.h"
#include "./core_shell.h"


namespace py = pybind11;

inline ScattererProperties make_scatterer_properties_from_material(
    const BaseSourceSet& source,
    py::object materials
)
{
    py::object BaseMaterial =
        py::module::import("PyOptik.material.base_class").attr("BaseMaterial");

    if (!py::isinstance<py::list>(materials) &&
        !py::isinstance<py::tuple>(materials))
        throw py::type_error("material must be a list of BaseMaterial");

    py::list items = py::cast<py::list>(materials);

    std::vector<std::vector<complex128>> spectral_values;
    std::vector<std::string> material_names;

    py::object wavelength = py::cast(source.wavelength);

    for (auto item : items)
    {
        if (!py::isinstance(item, BaseMaterial))
            throw py::type_error("All elements of material must be BaseMaterial");

        py::object result =
            item.attr("compute_refractive_index")(wavelength);

        spectral_values.push_back(
            cast_scalar_or_array_to_vector_complex128(result)
        );

        std::string name =
            py::cast<std::string>(item.attr("__repr__")());

        material_names.push_back(name);
    }

    return ScattererProperties(spectral_values, material_names);
}

inline MediumProperties make_medium_properties_from_material(
    const BaseSourceSet& source,
    py::object materials
)
{
    py::object BaseMaterial =
        py::module::import("PyOptik.material.base_class").attr("BaseMaterial");

    if (!py::isinstance<py::list>(materials) &&
        !py::isinstance<py::tuple>(materials))
        throw py::type_error("medium_material must be a list of BaseMaterial");

    py::list items = py::cast<py::list>(materials);

    std::vector<std::vector<double>> spectral_values;
    std::vector<std::string> material_names;

    py::object wavelength = py::cast(source.wavelength);

    for (auto item : items)
    {
        if (!py::isinstance(item, BaseMaterial))
            throw py::type_error("All elements of medium_material must be BaseMaterial");

        py::object result =
            item.attr("compute_refractive_index")(wavelength);

        spectral_values.push_back(
            cast_scalar_or_array_to_vector_double(result)
        );

        std::string name =
            py::cast<std::string>(item.attr("__repr__")());

        material_names.push_back(name);
    }

    return MediumProperties(spectral_values, material_names);
}




void register_scatterer_set(py::module& module) {
    py::object ureg = get_shared_ureg();


    py::class_<ScattererSet, std::shared_ptr<ScattererSet>>(module, "BaseScattererSet")
    .def_readonly("is_sequential", &ScattererSet::is_sequential)
    ;

    py::class_<SphereSet, ScattererSet, std::shared_ptr<SphereSet>>(module, "SphereSet")
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& refractive_index,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );
                    std::vector<complex128> refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            refractive_index.attr("to")("RIU").attr("magnitude")
                        );
                    std::vector<double> medium_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<SphereSet>(
                        diameter_value,
                        ScattererProperties(refractive_index_value),
                        MediumProperties(medium_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents an infinite cylinder scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding infinite cylinder scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding infinite cylinder-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the infinite cylinder scatterers in meters.
                refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the infinite cylinder scatterers themselves.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& refractive_index,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    std::vector<double> medium_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<SphereSet>(
                        diameter_value,
                        ScattererProperties(refractive_index_value),
                        MediumProperties(medium_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their refractive_index to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the spherical scatterers in meters.
                refractive_index : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the spherical scatterers themselves.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& material,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> medium_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<SphereSet>(
                        diameter_value,
                        make_scatterer_properties_from_material(source, material),
                        MediumProperties(medium_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the spherical scatterers in meters.
                material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the spherical scatterers themselves.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& refractive_index,
                    const py::object& medium_material
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<SphereSet>(
                        diameter_value,
                        ScattererProperties(refractive_index_value),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the spherical scatterers in meters.
                refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the spherical scatterers themselves.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& material,
                    const py::object& medium_material
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    return std::make_shared<SphereSet>(
                        diameter_value,
                        make_scatterer_properties_from_material(source, material),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the spherical scatterers in meters.
                material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the spherical scatterers themselves.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def_readonly_static("attributes", &SphereSet::attributes)
        .def_readonly_static("available_measure_list", &SphereSet::available_measure_list)
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
        .def_static(
            "build_sequential",
            [ureg](
                const size_t& target_size,
                const py::object& diameter,
                const py::object& refractive_index,
                const py::object& medium_refractive_index
            ) {
                std::vector<double> diameter_value =
                    cast_scalar_or_array_to_vector_double(diameter.attr("to")("meter").attr("magnitude"));

                std::vector<complex128> refractive_index_value =
                    cast_scalar_or_array_to_vector_complex128(refractive_index.attr("to")("RIU").attr("magnitude"));

                std::vector<double> medium_refractive_index_value =
                    cast_scalar_or_array_to_vector_double(medium_refractive_index.attr("to")("RIU").attr("magnitude"));

                diameter_value = broadcast_vector_double("diameter", diameter_value, target_size);
                refractive_index_value = broadcast_vector_complex128("refractive_index", refractive_index_value, target_size);
                medium_refractive_index_value = broadcast_vector_double("medium_refractive_index", medium_refractive_index_value, target_size);

                ScattererProperties scatterer_properties(refractive_index_value);
                MediumProperties medium_properties(medium_refractive_index_value);

                return std::make_shared<SphereSet>(
                    diameter_value,
                    scatterer_properties,
                    medium_properties,
                    true
                );
            },
            py::arg("total_size"),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                Construct a sequential SphereSet by broadcasting scalar or single element inputs.

                Parameters
                ----------
                total_size : int or None
                    Target size for broadcasting. If None, uses the size of the diameter vector after conversion.
                diameter : Quantity or array-like Quantity
                    Diameter(s) of the spheres. Scalars or length 1 arrays broadcast to total_size.
                refractive_index : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the spherical scatterers themselves.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

                Returns
                -------
                SphereSet
                    Instance with is_sequential = True.
            )pdoc"
        )
        ;
















    // Binding for CYLINDER::Set
    py::class_<InfiniteCylinderSet, ScattererSet, std::shared_ptr<InfiniteCylinderSet>>(module, "InfiniteCylinderSet")
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& refractive_index,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );
                    std::vector<complex128> refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            refractive_index.attr("to")("RIU").attr("magnitude")
                        );
                    std::vector<double> medium_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<InfiniteCylinderSet>(
                        diameter_value,
                        ScattererProperties(refractive_index_value),
                        MediumProperties(medium_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents an infinite cylinder scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding infinite cylinder scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding infinite cylinder-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the infinite cylinder scatterers in meters.
                refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the infinite cylinder scatterers themselves.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& material,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> medium_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<InfiniteCylinderSet>(
                        diameter_value,
                        make_scatterer_properties_from_material(source, material),
                        MediumProperties(medium_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents an infinite cylinder scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding infinite cylinder scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding infinite cylinder-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the infinite cylinder scatterers in meters.
                material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the infinite cylinder scatterers themselves.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& refractive_index,
                    const py::object& medium_material
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<InfiniteCylinderSet>(
                        diameter_value,
                        ScattererProperties(refractive_index_value),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents an infinite cylinder scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding infinite cylinder scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding infinite cylinder-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the infinite cylinder scatterers in meters.
                refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the infinite cylinder scatterers themselves.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& diameter,
                    const py::object& material,
                    const py::object& medium_material
                ) {

                    std::vector<double> diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            diameter.attr("to")("meter").attr("magnitude")
                        );

                    return std::make_shared<InfiniteCylinderSet>(
                        diameter_value,
                        make_scatterer_properties_from_material(source, material),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents an infinite cylinder scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding infinite cylinder scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding infinite cylinder-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                diameter : Length
                    Diameter(s) of the infinite cylinder scatterers in meters.
                material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the infinite cylinder scatterers themselves.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def_readonly_static("attributes", &InfiniteCylinderSet::attributes)
        .def_readonly_static("available_measure_list", &InfiniteCylinderSet::available_measure_list)
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
        .def_static(
            "build_sequential",
            [ureg](
                const size_t& target_size,
                const py::object& diameter,
                const py::object& refractive_index,
                const py::object& medium_refractive_index
            ) {
                std::vector<double> diameter_value =
                    cast_scalar_or_array_to_vector_double(diameter.attr("to")("meter").attr("magnitude"));

                std::vector<complex128> refractive_index_value =
                    cast_scalar_or_array_to_vector_complex128(refractive_index.attr("to")("RIU").attr("magnitude"));

                std::vector<double> medium_refractive_index_value =
                    cast_scalar_or_array_to_vector_double(medium_refractive_index.attr("to")("RIU").attr("magnitude"));

                diameter_value = broadcast_vector_double("diameter", diameter_value, target_size);
                refractive_index_value = broadcast_vector_complex128("refractive_index", refractive_index_value, target_size);
                medium_refractive_index_value = broadcast_vector_double("medium_refractive_index", medium_refractive_index_value, target_size);

                ScattererProperties scatterer_properties(refractive_index_value);
                MediumProperties medium_properties(medium_refractive_index_value);

                return std::make_shared<InfiniteCylinderSet>(
                    diameter_value,
                    scatterer_properties,
                    medium_properties,
                    true
                );
            },
            py::arg("total_size"),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                Construct a sequential InfiniteCylinderSet by broadcasting scalar or single element inputs.

                Parameters
                ----------
                total_size : int or None
                    Target size for broadcasting. If None, uses the size of the diameter vector after conversion.
                diameter : Quantity or array-like Quantity
                    Diameter(s) of the cylinders. Scalars or length 1 arrays broadcast to total_size.
                refractive_index : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the cylindrical scatterers themselves.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

                Returns
                -------
                InfiniteCylinderSet
                    Instance with is_sequential = True.
            )pdoc"
        )
        ;

















    // Binding for CORESHELL::Set
    py::class_<CoreShellSet, ScattererSet, std::shared_ptr<CoreShellSet>>(module, "CoreShellSet")
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_refractive_index,
                    const py::object& shell_refractive_index,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> core_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            core_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    std::vector<complex128> shell_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            shell_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    std::vector<double> medium_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        ScattererProperties(core_refractive_index_value),
                        ScattererProperties(shell_refractive_index_value),
                        MediumProperties(medium_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents a core-shell scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the core of the spherical scatterers.
                shell_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the shell of the spherical scatterers.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_refractive_index,
                    const py::object& shell_refractive_index,
                    const py::object& medium_material
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> core_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            core_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    std::vector<complex128> shell_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            shell_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        ScattererProperties(core_refractive_index_value),
                        ScattererProperties(shell_refractive_index_value),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_refractive_index"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the core of the spherical scatterers.
                shell_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the shell of the spherical scatterers.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_refractive_index,
                    const py::object& shell_material,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> core_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            core_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    std::vector<double> medium_refractive_index_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        ScattererProperties(core_refractive_index_value),
                        make_scatterer_properties_from_material(source, shell_material),
                        MediumProperties(medium_refractive_index_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_material"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents a core-shell scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding core-shell scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding core-shell-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the core of the spherical scatterers.
                shell_material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the shell of the spherical scatterers.
                medium_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the medium surrounding the scatterers.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_refractive_index,
                    const py::object& shell_material,
                    const py::object& medium_material
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> core_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            core_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        ScattererProperties(core_refractive_index_value),
                        make_scatterer_properties_from_material(source, shell_material),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_material"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents a core-shell scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding core-shell scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding core-shell-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the core of the spherical scatterers.
                shell_material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the shell of the spherical scatterers.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_material,
                    const py::object& shell_refractive_index,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> shell_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            shell_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    std::vector<double> medium_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        make_scatterer_properties_from_material(source, core_material),
                        ScattererProperties(shell_refractive_index_value),
                        MediumProperties(medium_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_material"),
            py::arg("shell_refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents a core-shell scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_material : List[BaseMaterial]]
                    Refractive index or indices of the core of the spherical scatterers.
                shell_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the shell of the spherical scatterers.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_material,
                    const py::object& shell_refractive_index,
                    const py::object& medium_material
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<complex128> shell_refractive_index_value =
                        cast_scalar_or_array_to_vector_complex128(
                            shell_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        make_scatterer_properties_from_material(source, core_material),
                        ScattererProperties(shell_refractive_index_value),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_material"),
            py::arg("shell_refractive_index"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents a spherical scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding spherical scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding spherical-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_material : List[BaseMaterial]
                    Refractive index or indices of the core of the spherical scatterers.
                shell_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the shell of the spherical scatterers.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_material,
                    const py::object& shell_material,
                    const py::object& medium_refractive_index
                ) {

                    std::vector<double> core_diameter_value =
                        cast_scalar_or_array_to_vector_double(
                            core_diameter.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> shell_thickness_value =
                        cast_scalar_or_array_to_vector_double(
                            shell_thickness.attr("to")("meter").attr("magnitude")
                        );

                    std::vector<double> medium_refractive_index_value =
                        cast_scalar_or_array_to_vector_double(
                            medium_refractive_index.attr("to")("RIU").attr("magnitude")
                        );

                    return std::make_shared<CoreShellSet>(
                        core_diameter_value,
                        shell_thickness_value,
                        make_scatterer_properties_from_material(source, core_material),
                        make_scatterer_properties_from_material(source, shell_material),
                        MediumProperties(medium_refractive_index_value),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_material"),
            py::arg("shell_material"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                A data class that represents a core-shell scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding core-shell scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding core-shell-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_material : List[BaseMaterial]
                    Material or materials of the core of the spherical scatterers.
                shell_material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the shell of the spherical scatterers.
                medium_refractive_index : List[RefractiveIndex]
                    Refractive index or indices of the medium surrounding the scatterers.
            )pdoc"
        )
        .def(
            py::init(
                [ureg](
                    const BaseSourceSet& source,
                    const py::object& core_diameter,
                    const py::object& shell_thickness,
                    const py::object& core_material,
                    const py::object& shell_material,
                    const py::object& medium_material
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
                        make_scatterer_properties_from_material(source, core_material),
                        make_scatterer_properties_from_material(source, shell_material),
                        make_medium_properties_from_material(source, medium_material),
                        false
                    );
                }
            ),
            py::arg("source"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_material"),
            py::arg("shell_material"),
            py::arg("medium_material"),
            R"pdoc(
                A data class that represents a core-shell scatterer configuration used in PyMieSim simulations.

                This class provides specific implementations for setting up and binding core-shell scatterers
                with their material to a simulation environment. It extends the `BaseScatterer` class by
                adding core-shell-specific attributes and methods for handling simulation setups.

                Parameters
                ----------
                source : PyMieSim.binary.interface_experiment.BaseSourceSet
                    Light source configuration for the simulation.
                core_diameter : Length
                    Diameter(s) of the core of the spherical scatterers in meters.
                shell_thickness : Length
                    Thickness of the shell of the spherical scatterers in meters.
                core_material : List[BaseMaterial]
                    Material or materials of the core of the spherical scatterers.
                shell_material : List[BaseMaterial] | List[RefractiveIndex]
                    Material or materials of the shell of the spherical scatterers.
                medium_material : List[BaseMaterial]
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.
            )pdoc"
        )
        .def_readonly_static("attributes", &CoreShellSet::attributes)
        .def_readonly_static("available_measure_list", &CoreShellSet::available_measure_list)
        .def(
            "get_mapping",
            [ureg](const CoreShellSet& self) {
                py::dict mapping;
                mapping["scatterer:core_diameter"] = py::cast(self.core_diameter) * ureg.attr("meter");
                mapping["scatterer:shell_thickness"] = py::cast(self.shell_thickness) * ureg.attr("meter");
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
        .def_static(
            "build_sequential",
            [ureg](
                const size_t& target_size,
                const py::object& core_diameter,
                const py::object& shell_thickness,
                const py::object& core_refractive_index,
                const py::object& shell_refractive_index,
                const py::object& medium_refractive_index
            ) {
                std::vector<double> core_diameter_value = \
                    cast_scalar_or_array_to_vector_double(core_diameter.attr("to")("meter").attr("magnitude"));

                std::vector<double> shell_thickness_value = \
                    cast_scalar_or_array_to_vector_double(shell_thickness.attr("to")("meter").attr("magnitude"));

                std::vector<complex128> core_refractive_index_value =
                    cast_scalar_or_array_to_vector_complex128(core_refractive_index.attr("to")("RIU").attr("magnitude"));

                std::vector<complex128> shell_refractive_index_value =
                    cast_scalar_or_array_to_vector_complex128(shell_refractive_index.attr("to")("RIU").attr("magnitude"));

                std::vector<double> medium_refractive_index_value =
                    cast_scalar_or_array_to_vector_double(medium_refractive_index.attr("to")("RIU").attr("magnitude"));

                core_diameter_value = broadcast_vector_double("core_diameter", core_diameter_value, target_size);
                shell_thickness_value = broadcast_vector_double("shell_thickness", shell_thickness_value, target_size);
                core_refractive_index_value = broadcast_vector_complex128("core_refractive_index", core_refractive_index_value, target_size);
                shell_refractive_index_value = broadcast_vector_complex128("shell_refractive_index", shell_refractive_index_value, target_size);
                medium_refractive_index_value = broadcast_vector_double("medium_refractive_index", medium_refractive_index_value, target_size);

                return std::make_shared<CoreShellSet>(
                    core_diameter_value,
                    shell_thickness_value,
                    core_refractive_index_value,
                    shell_refractive_index_value,
                    medium_refractive_index_value,
                    true
                );
            },
            py::arg("total_size"),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_refractive_index"),
            py::arg("medium_refractive_index"),
            R"pdoc(
                Construct a sequential CoreShellSet by broadcasting scalar or single element inputs.

                Parameters
                ----------
                total_size : int or None
                    Target size for broadcasting. If None, uses the size of the core_diameter vector after conversion.
                core_diameter : Quantity or array-like Quantity
                    Diameter(s) of the core. Scalars or length 1 arrays broadcast to total_size.
                shell_thickness : Quantity or array-like Quantity
                    Thickness(es) of the shell. Scalars or length 1 arrays broadcast to total_size.
                core_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the core.
                shell_refractive_index : List[BaseMaterial] | List[RefractiveIndex]
                    Refractive index or indices of the shell.
                medium_refractive_index : List, optional
                    BaseMaterial(s) defining the medium, used if `medium_index` is not provided.

                Returns
                -------
                CoreShellSet
                    Instance with is_sequential = True.
            )pdoc"
        )
        ;
}

// -
