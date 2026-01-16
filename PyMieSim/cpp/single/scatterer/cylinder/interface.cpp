#include <pybind11/pybind11.h>
#include "./cylinder.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>

namespace py = pybind11;

void register_cylinder(py::module_& module) {

    // Binding for Cylinder class
    py::class_<Cylinder, BaseScatterer, std::shared_ptr<Cylinder>>(module, "CYLINDER")
        .def(
            "__init__",
            [](Cylinder &instance, py::object diameter, py::object refractive_index, py::object medium_refractive_index, std::shared_ptr<BaseSource> source) {
                py::object ureg = get_shared_ureg();

                double diameter_meter = (diameter.attr("to")(ureg.attr("meter"))).attr("magnitude").cast<double>();
                std::complex<double> refractive_index_value = refractive_index.attr("to")(ureg.attr("RIU")).attr("magnitude").cast<std::complex<double>>();
                double medium_refractive_index_value = (medium_refractive_index.attr("to")(ureg.attr("RIU"))).attr("magnitude").cast<double>();
                new (&instance) Cylinder(diameter_meter, refractive_index_value, medium_refractive_index_value, std::move(source));
            },
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            R"pbdoc(
                Constructor for CYLINDER, initializing it with physical and optical properties.

                Parameters
                ----------
                diameter : float
                    The diameter of the cylinder.
                refractive_index : complex
                    The refractive index of the cylinder.
                medium_refractive_index : float
                    The refractive index of the surrounding medium.
                source : BaseSource
                    The source of the incident light.
            )pbdoc"
        )
        .def_readonly(
            "source",
            &Cylinder::source,
            "Source of the cylinder."
        )
        .def_property_readonly(
            "diameter",
            [&](const Cylinder &self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.diameter) * ureg.attr("meter")).attr("to_compact")();
            },
            "Diameter of the cylinder."
        )
        .def_property_readonly(
            "refractive_index",
            [](const Cylinder &self) {
                py::object ureg = get_shared_ureg();

                py::object magnitude = py::cast(self.refractive_index);

                return (magnitude * ureg.attr("RIU")).attr("to_compact")();
            },
            "Refractive index of the cylinder."
        )
        .def_property_readonly(
            "medium_refractive_index",
            [](const Cylinder &self) {
                py::object ureg = get_shared_ureg();
                py::object magnitude = py::cast(self.medium_refractive_index);

                return (magnitude * ureg.attr("RIU")).attr("to_compact")();
            },
            "Refractive index of the surrounding medium."
        )
        .def_property_readonly(
            "radius",
            [&](const Cylinder &self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.diameter / 2.0) * ureg.attr("meter")).attr("to_compact")();
            },
            "Radius of the cylinder."
        )
        .def_property("a1n",
            [](Cylinder& self) {return vector_as_numpy_view(self, self.a1n);},
            [](Cylinder& self,
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
            [](Cylinder& self) {return vector_as_numpy_view(self, self.b1n);},
            [](Cylinder& self,
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
            [](Cylinder& self) {return vector_as_numpy_view(self, self.a2n);},
            [](Cylinder& self,
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
            [](Cylinder& self) {return vector_as_numpy_view(self, self.b2n);},
            [](Cylinder& self,
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
