#include <pybind11/pybind11.h>
#include "sphere.h"


void register_sphere(pybind11::module_& module) {

    pybind11::class_<Sphere, BaseScatterer>(module, "SPHERE")
        .def(
            pybind11::init<const double, const complex128, const double, const BaseSource&, size_t>(),
            pybind11::arg("diameter"),
            pybind11::arg("refractive_index"),
            pybind11::arg("medium_refractive_index"),
            pybind11::arg("source"),
            pybind11::arg("max_order") = 0,
            R"pbdoc(
                Constructor for SPHERE, initializing it with physical and optical properties.

                Parameters
                ----------
                diameter : float
                    The diameter of the sphere.
                refractive_index : complex
                    The refractive index of the sphere.
                medium_refractive_index : float
                    The refractive index of the surrounding medium.
                source : BaseSource
                    The source of the incident light.
                max_order : int, optional
                    The maximum order of spherical harmonics to use in the scattering calculation. Default is 0.
            )pbdoc"
        )
        .def_property("an",
            [](Sphere& self) {return vector_as_numpy_view(self, self.an);},
            [](Sphere& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
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
            [](Sphere& self) {return vector_as_numpy_view(self, self.bn);},
            [](Sphere& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
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
            [](Sphere& self) {return vector_as_numpy_view(self, self.cn);},
            [](Sphere& self,
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
            [](Sphere& self) {return vector_as_numpy_view(self, self.dn);},
            [](Sphere& self,
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
