#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include "./sphere.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>

namespace py = pybind11;

void register_sphere(py::module_& module) {
    py::object ureg = get_shared_ureg();

    py::class_<Sphere, BaseScatterer, std::shared_ptr<Sphere>>(module, "Sphere")
        .def(
            py::init([ureg](
                py::object diameter,
                py::object refractive_index,
                py::object medium_refractive_index,
                std::shared_ptr<BaseSource> source,
                std::size_t max_order
            ) {
                py::object units_length = py::module_::import("PyMieSim.units").attr("Length");
                py::object units_riu = py::module_::import("PyMieSim.units").attr("RefractiveIndex");

                diameter = units_length.attr("check")(diameter);
                refractive_index = units_riu.attr("check")(refractive_index);
                medium_refractive_index = units_riu.attr("check")(medium_refractive_index);

                double diameter_meter = diameter.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                complex128 refractive_index_value = refractive_index.attr("to")(ureg.attr("refractive_index_unit")).attr("magnitude").cast<complex128>();

                double medium_refractive_index_value = medium_refractive_index.attr("to")(ureg.attr("refractive_index_unit")).attr("magnitude").cast<double>();

                return std::make_shared<Sphere>(
                    diameter_meter,
                    refractive_index_value,
                    medium_refractive_index_value,
                    std::move(source),
                    max_order
                );
            }),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Constructor for SPHERE, initializing it with physical and optical properties.

                Parameters
                ----------
                diameter : float[Length]
                    The diameter of the sphere.
                refractive_index : complex[RefractiveIndex]
                    The refractive index of the sphere.
                medium_refractive_index : float[RefractiveIndex]
                    The refractive index of the surrounding medium.
                source : BaseSource
                    The source of the incident light.
                max_order : int, optional
                    The maximum order of spherical harmonics to use in the scattering calculation. Default is 0.
            )pbdoc"
        )
        .def_readonly(
            "source",
            &Sphere::source,
            "Source of the sphere."
        )
        .def_readonly(
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
        .def_property_readonly(
            "refractive_index",
            [ureg](const Sphere &self) {
                py::object magnitude = py::cast(self.refractive_index);

                return (magnitude * ureg.attr("RIU")).attr("to_compact")();
            },
            "Refractive index of the sphere."
        )
        .def_property_readonly(
            "medium_refractive_index",
            [ureg](const Sphere &self) {
                py::object magnitude = py::float_(self.medium_refractive_index);

                return (magnitude * ureg.attr("RIU")).attr("to_compact")();
            },
            "Refractive index of the surrounding medium."
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
        .def_property("an",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.an);},
            [ureg](Sphere& self, py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr)
            {
                vector_assign_from_numpy(self.an, arr);
            },
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
        .def_property("bn",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.bn);},
            [ureg](Sphere& self, py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr)
            {
                vector_assign_from_numpy(self.bn, arr);
            },
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
        .def_property("cn",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.cn);},
            [ureg](Sphere& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.cn, arr);
            },
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
        .def_property("dn",
            [ureg](Sphere& self) {return vector_as_numpy_view(self, self.dn);},
            [ureg](Sphere& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.dn, arr);
            },
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
            )pbdoc")
        ;
}
