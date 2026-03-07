#include <pybind11/pybind11.h>
#include "./cylinder.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>
#include <single/scatterer/utils.h>

namespace py = pybind11;

void register_cylinder(py::module_& module) {
    py::object ureg = get_shared_ureg();

    // Binding for InfiniteCylinder class
    py::class_<InfiniteCylinder, BaseScatterer, std::shared_ptr<InfiniteCylinder>>(module, "InfiniteCylinder")
        .def(
            py::init([ureg](
                py::object diameter,
                const py::object& material,
                const py::object& medium,
                const std::shared_ptr<BaseSource> source,
                std::size_t max_order
            ) {
                double diameter_meter = diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                const std::shared_ptr<BaseMaterial> parsed_material =
                    parse_material_object(material, ureg);

                const std::shared_ptr<BaseMedium> parsed_medium =
                    parse_medium_object(medium, ureg);

                return std::make_shared<InfiniteCylinder>(
                    diameter_meter,
                    std::move(parsed_material),
                    std::move(parsed_medium),
                    std::move(source),
                    max_order
                );
            }),
            py::arg("diameter"),
            py::arg("material"),
            py::arg("medium"),
            py::arg("source"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Constructor for CYLINDER, initializing it with physical and optical properties.

                Parameters
                ----------
                diameter : float
                    The diameter of the cylinder.
                material : BaseMaterial
                    The material of the cylinder.
                medium : BaseMedium
                    The surrounding medium of the cylinder.
                source : BaseSource
                    The source of the incident light.
            )pbdoc"
        )
        .def_readonly(
            "source",
            &InfiniteCylinder::source,
            "Source of the cylinder."
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
        .def_property("a1n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.a1n);},
            [ureg](InfiniteCylinder& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.a1n, arr);
            },
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
        .def_property("b1n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.b1n);},
            [ureg](InfiniteCylinder& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.b1n, arr);
            },
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
        .def_property("a2n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.a2n);},
            [ureg](InfiniteCylinder& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.a2n, arr);
            },
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
        .def_property("b2n",
            [ureg](InfiniteCylinder& self) {return vector_as_numpy_view(self, self.b2n);},
            [ureg](InfiniteCylinder& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.b2n, arr);
            },
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
    ;

}
