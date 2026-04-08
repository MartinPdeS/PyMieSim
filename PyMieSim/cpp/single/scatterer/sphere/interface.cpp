#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include "./sphere.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>
#include <single/scatterer/utils.h>

namespace py = pybind11;


void register_sphere(py::module_& module) {
    py::object ureg = get_shared_ureg();

    py::class_<Sphere, BaseScatterer, std::shared_ptr<Sphere>>(
            module,
            "Sphere",
            R"pbdoc(
                A class representing a spherical scatterer, defined by its diameter, material, and surrounding medium.

                The Sphere class provides methods to compute the scattering coefficients based on Mie theory, as well as properties to access its physical and optical characteristics.
            )pbdoc"
        )
        .def(
            py::init([ureg](
                const py::object& diameter,
                const py::object& material,
                const py::object& medium,
                const std::size_t max_order
            ) {
                const double diameter_meter = diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                const std::shared_ptr<BaseMaterial> parsed_material =
                    parse_material_object(material, ureg);

                const std::shared_ptr<BaseMedium> parsed_medium =
                    parse_medium_object(medium, ureg);

                return std::make_shared<Sphere>(
                    diameter_meter,
                    std::move(parsed_material),
                    std::move(parsed_medium),
                    max_order
                );
            }),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Constructor for Sphere, initializing it with physical and optical properties.

                Parameters
                ----------
                diameter : float[Length]
                    The diameter of the sphere.
                material : BaseMaterial or complex[RefractiveIndex]
                    The material of the sphere, or a constant complex refractive index.
                medium : BaseMedium or float[RefractiveIndex]
                    The surrounding medium, or a constant real refractive index.
                max_order : int, optional
                    The maximum order of spherical harmonics to use in the scattering calculation. Default is 0.
            )pbdoc"
        )
        .def_readonly_static(
            "property_names",
            &Sphere::property_names,
            "Property names of the sphere."
        )
        .def(
            "print_properties",
            &Sphere::print_properties,
            "Prints the properties of the sphere."
        )
        .def_property_readonly(
            "diameter",
            [ureg](const Sphere &self) {
                return (py::float_(self.diameter) * ureg.attr("meter")).attr("to_compact")();
            },
            "Diameter of the sphere."
        )
        .def_readonly(
            "material",
            &Sphere::material,
            "Material of the sphere."
        )
        .def_readonly(
            "medium",
            &Sphere::medium,
            "Medium surrounding the sphere."
        )
        .def_property_readonly(
            "radius",
            [ureg](const Sphere &self) {

                return (py::float_(self.diameter / 2.0) * ureg.attr("meter")).attr("to_compact")();
            },
            "Radius of the sphere."
        )
        .def_property_readonly(
            "volume",
            [ureg](const Sphere &self) {
                double volume = (4.0 / 3.0) * Constants::PI * std::pow(self.diameter / 2.0, 3);
                return (py::float_(volume) * ureg.attr("meter**3")).attr("to_compact")();
            },
            "Volume of the sphere."
        )
        .def_property_readonly(
            "an",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.an);},
            R"pbdoc(
                Returns :math:`a_n` coefficient as defined in Eq:III.88 of B&B:

                .. math::
                    a_n = \frac{
                        \mu_{sp} \Psi_n(\alpha) \Psi_n^\prime(\beta) - \mu M \Psi_n^\prime(\alpha) \Psi_n(\beta)
                    }{
                        \mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta) - \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)
                    }

                With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)


                Returns
                -------
                list
                    A list of 'an' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def_property_readonly(
            "bn",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.bn);},
            R"pbdoc(
                Returns :math:`b_n` coefficient as defined in Eq:III.89 of B&B:

                .. math::
                    b_n = \frac{
                        \mu M \Psi_n(\alpha) \Psi_n^\prime(\beta) - \mu_{sp} \Psi_n^\prime(\alpha) \Psi_n(\beta)
                    }{
                        \mu M \xi_n(\alpha) \Psi_n^\prime(\beta) - \mu_{sp} \xi_n^\prime (\alpha) \Psi_n(\beta)
                    }

                With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

                Returns
                -------
                list
                    A list of 'bn' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def_property_readonly(
            "cn",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.cn);},
            R"pbdoc(
                For future purpose only!
                Returns :math:`c_n` coefficient as defined in Eq:III.90 of B&B:

                .. math::
                    c_n = \frac{
                        \mu_{sp} M \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) - \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]
                    }{
                        \mu_{sp} \xi_n(\alpha) \Psi_n^\prime(\beta) - \mu M \xi_n^\prime (\alpha) \Psi_n(\beta)
                    }

                With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

                Returns
                -------
                list
                    A list of 'cn' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def_property_readonly(
            "dn",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.dn);},
            R"pbdoc(
                For future purpose only!
                Returns :math:`d_n` coefficient as defined in Eq:III.91 of B&B:

                .. math::
                    d_n = \frac{
                        \mu M^2 \big[ \xi_n(\alpha) \Psi_n^\prime(\alpha) - \xi_n^\prime(\alpha) \Psi_n(\alpha) \big]
                    }{
                        \mu M \xi_n(\alpha) \Psi_n^\prime(\beta) - \mu_{sp} M \xi_n^\prime (\alpha) \Psi_n(\beta)
                    }

                With :math:`M = \frac{k_{sp}}{k}` (Eq:I.103)

                Returns
                -------
                list
                    A list of 'dn' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        .def(
            "add_to_scene",
            [](const Sphere& self, const py::object& scene, py::kwargs kwargs) {
                py::module_ pyvista = py::module_::import("pyvista");

                py::object physical_sphere = pyvista.attr("Sphere")(
                    py::arg("radius") = py::float_(0.1),
                    py::arg("center") = py::make_tuple(0.0, 0.0, 0.0),
                    py::arg("theta_resolution") = 50,
                    py::arg("phi_resolution") = 50
                );

                py::object unit_sphere = pyvista.attr("Sphere")(
                    py::arg("radius") = py::float_(1.0),
                    py::arg("center") = py::make_tuple(0.0, 0.0, 0.0),
                    py::arg("theta_resolution") = 50,
                    py::arg("phi_resolution") = 50
                );

                scene.attr("add_mesh")(physical_sphere, **kwargs);

                py::dict unit_sphere_kwargs;
                unit_sphere_kwargs["opacity"] = py::float_(0.15);

                if (kwargs.contains("color")) {
                    unit_sphere_kwargs["color"] = kwargs["color"];
                }
                else {
                    unit_sphere_kwargs["color"] = py::str("white");
                }

                scene.attr("add_mesh")(unit_sphere, **unit_sphere_kwargs);
            },
            py::arg("scene"),
            R"pbdoc(
                Adds the sphere to a given PyVista scene for visualization.

                This method creates a physical representation of the sphere based on its diameter and adds it to the provided PyVista scene. It also adds a unit sphere for reference, which can be customized in appearance using keyword arguments.

                Parameters
                ----------
                scene : pyvista.Plotter
                    The PyVista Plotter object to which the sphere will be added.
                **kwargs
                    Additional keyword arguments to customize the appearance of the physical sphere and the unit sphere. For example, you can specify 'color' to set the color of both spheres, and 'opacity' to set their transparency.
            )pbdoc"
        )
        ;
}
