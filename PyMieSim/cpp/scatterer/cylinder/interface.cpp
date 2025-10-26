#include <pybind11/pybind11.h>
#include "cylinder.h"
#include "utils/numpy_interface.h"



void register_cylinder(pybind11::module_& module) {

    // Binding for Cylinder class
    pybind11::class_<Cylinder, BaseScatterer>(module, "CYLINDER")
        .def(
            pybind11::init<double, complex128, double, BaseSource&>(),
            pybind11::arg("diameter"),
            pybind11::arg("refractive_index"),
            pybind11::arg("medium_refractive_index"),
            pybind11::arg("source"),
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
        .def_property("a1n",
            [](Cylinder& self) {return vector_as_numpy_view(self, self.a1n);},
            [](Cylinder& self,
               pybind11::array_t<std::complex<double>, pybind11::array::c_style | pybind11::array::forcecast> arr) {
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
               pybind11::array_t<std::complex<double>, pybind11::array::c_style | pybind11::array::forcecast> arr) {
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
               pybind11::array_t<std::complex<double>, pybind11::array::c_style | pybind11::array::forcecast> arr) {
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
               pybind11::array_t<std::complex<double>, pybind11::array::c_style | pybind11::array::forcecast> arr) {
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
