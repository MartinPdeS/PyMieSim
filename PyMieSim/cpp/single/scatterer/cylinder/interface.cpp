#include <pybind11/pybind11.h>
#include "./cylinder.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>
#include <single/scatterer/utils.h>

namespace py = pybind11;

void register_cylinder(py::module_& module) {
    py::object ureg = get_shared_ureg();

    // Binding for InfiniteCylinder class
    py::class_<InfiniteCylinder, BaseScatterer, std::shared_ptr<InfiniteCylinder>>(
            module,
            "InfiniteCylinder",
            R"pbdoc(
                A class representing an infinite cylindrical scatterer, defined by its diameter, material, and surrounding medium.

                The InfiniteCylinder class provides methods to compute the scattering coefficients based on Mie theory for cylinders, as well as properties to access its physical and optical characteristics.
            )pbdoc"
        )
        .def(
            py::init(
                [ureg](
                    py::object diameter,
                    const py::object& material,
                    const py::object& medium,
                    std::size_t max_order
                ) {
                    double diameter_meter = diameter.attr("to")("meter").attr("magnitude").cast<double>();

                    const std::shared_ptr<BaseMaterial> parsed_material =
                        parse_material_object(material, ureg);

                    const std::shared_ptr<BaseMedium> parsed_medium =
                        parse_medium_object(medium, ureg);

                    return std::make_shared<InfiniteCylinder>(
                        diameter_meter,
                        std::move(parsed_material),
                        std::move(parsed_medium),
                        max_order
                    );
                }
            ),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Constructor for InfiniteCylinder, initializing it with physical and optical properties.

                Parameters
                ----------
                diameter : float
                    The diameter of the cylinder.
                material : BaseMaterial
                    The material of the cylinder.
                medium : BaseMedium
                    The surrounding medium of the cylinder.
                max_order : int, optional
                    The maximum order of the scattering coefficients to compute (default is 0, which means it will be determined automatically based on the size parameter).
            )pbdoc"
        )
        .def_readonly_static(
            "property_names",
            &InfiniteCylinder::property_names,
            "Property names of the cylinder."
        )
        .def_property_readonly(
            "diameter",
            [&](const InfiniteCylinder &self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.diameter) * ureg.attr("meter")).attr("to_compact")();
            },
            "Diameter of the cylinder."
        )
        .def(
            "print_properties",
            &InfiniteCylinder::print_properties,
            "Prints the properties of the cylinder."
        )
        .def_readonly(
            "material",
            &InfiniteCylinder::material,
            "Material of the cylinder."
        )
        .def_readonly(
            "medium",
            &InfiniteCylinder::medium,
            "Medium of the cylinder."
        )
        .def_property_readonly(
            "radius",
            [&](const InfiniteCylinder &self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.diameter / 2.0) * ureg.attr("meter")).attr("to_compact")();
            },
            "Radius of the cylinder."
        )
        .def_property_readonly(
            "a1n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.a1n);},
            R"pbdoc(
                Returns :math:`a_{1n}` coefficient as defined in ref[5]:

                .. math::
                    a_{1n} = \frac{
                        m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x)
                    }{
                        m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x)
                    }

                With :math:`m` being the refractive index of the medium and
                :math:`m_t` being the refractive index of the scatterer.

                Returns
                -------
                list
                    A list of a1n scattering coefficients for the cylinder. These coefficients describe the scattering behavior of the cylinder in the first mode of the spherical wave expansion.
            )pbdoc"
        )
        .def_property_readonly(
            "b1n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.b1n);},
            R"pbdoc(
                Returns :math:`b_{1n}` coefficient as defined in ref[5]:

                .. math::
                    b_{1n} = \frac{
                        m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x)
                    }{
                        m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x)
                    }

                With :math:`m` being the refractive index of the medium and
                :math:`m_t` being the refractive index of the scatterer.

                Returns
                -------
                list
                    A list of b1n scattering coefficients for the cylinder. These coefficients describe the scattering behavior of the cylinder in the first mode of the spherical wave expansion.
            )pbdoc"
        )
        .def_property_readonly(
            "a2n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.a2n);},
            R"pbdoc(
                Returns :math:`a_{2n}` coefficient as defined in ref[5]:

                .. math::
                    a_{2n} = \frac{
                        m_t J_n(m_t x) J_n^\prime (m x) - m J_n^\prime (m_t x) J_n(m x)
                    }{
                        m_t J_n(m_t x) H_n^\prime (m x) - m J_n^\prime (m_t x) H_n(m x)
                    }

                With :math:`m` being the refractive index of the medium and
                :math:`m_t` being the refractive index of the scatterer.

                Returns
                -------
                list
                    A list of a2n scattering coefficients for the cylinder. These coefficients describe the scattering behavior of the cylinder in the second mode of the spherical wave expansion.
            )pbdoc"
        )
        .def_property_readonly(
            "b2n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.b2n);},
            R"pbdoc(
                Returns :math:`b_{2n}` coefficient as defined in ref[5]:

                .. math::
                    b_{2n} = \frac{
                        m J_n(m_t x) J_n^\prime (m x) - m_t J_n^\prime (m_t x) J_n(m x)
                    }{
                        m J_n(m_t x) H_n^\prime (m x) - m_t J_n^\prime (m_t x) H_n(m x)
                    }

                With :math:`m` being the refractive index of the medium and
                :math:`m_t` being the refractive index of the scatterer.

                Returns
                -------
                list
                    A list of b2n scattering coefficients for the cylinder. These coefficients describe the scattering behavior of the cylinder in the second mode of the spherical wave expansion.
            )pbdoc"
        )
        .def(
            "add_to_scene",
            [](const InfiniteCylinder& self, const py::object& scene, py::kwargs kwargs) {
                py::module_ pyvista = py::module_::import("pyvista");

                py::object physical_cylinder = pyvista.attr("Cylinder")(
                    py::arg("center") = py::make_tuple(0.0, 0.0, 0.0),
                    py::arg("direction") = py::make_tuple(0.0, -1.0, 0.0),
                    py::arg("radius") = py::float_(0.1),
                    py::arg("height") = py::float_(2.0),
                    py::arg("resolution") = 100
                );

                py::object unit_sphere = pyvista.attr("Sphere")(
                    py::arg("radius") = py::float_(1.0),
                    py::arg("center") = py::make_tuple(0.0, 0.0, 0.0),
                    py::arg("theta_resolution") = 50,
                    py::arg("phi_resolution") = 50
                );

                scene.attr("add_mesh")(physical_cylinder, **kwargs);

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
                Adds a visual representation of the infinite cylinder to a given PyVista scene.

                This method creates a physical representation of the infinite cylinder using PyVista's Cylinder mesh and adds it to the provided scene. Additionally, it adds a semi-transparent unit sphere at the origin to help visualize the scale of the cylinder.

                Parameters
                ----------
                scene : pyvista.Plotter
                    The PyVista scene (Plotter) to which the cylinder will be added.
                **kwargs
                    Additional keyword arguments that can be passed to customize the appearance of the cylinder in the scene (e.g., color, opacity).
            )pbdoc"
        )
    ;

}
