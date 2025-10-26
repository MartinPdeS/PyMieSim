#include <pybind11/pybind11.h>
#include "coreshell.h"


void register_coreshell(pybind11::module_& module) {

    // Binding for CoreShell class
    pybind11::class_<CoreShell, BaseScatterer>(module, "CORESHELL")
        .def(pybind11::init<double, double, std::complex<double>, std::complex<double>, double, BaseSource&>(),
            pybind11::arg("core_diameter"),
            pybind11::arg("shell_thickness"),
            pybind11::arg("core_refractive_index"),
            pybind11::arg("shell_refractive_index"),
            pybind11::arg("medium_refractive_index"),
            pybind11::arg("source"),
            R"pbdoc(
                Constructor for CORESHELL, initializing it with physical and optical properties.

                Parameters
                ----------
                core_diameter : float
                    The diameter of the core of the spherical shell.
                shell_thickness : float
                    The thickness of the shell surrounding the core.
                core_refractive_index : complex
                    The refractive index of the core.
                shell_refractive_index : complex
                    The refractive index of the shell.
                medium_refractive_index : float
                    The refractive index of the surrounding medium.
                source : BaseSource
                    The source of the incident light.
            )pbdoc"
        )
        .def_property("an",
            [](Sphere& self) {return vector_as_numpy_view(self, self.an);},
            [](Sphere& self,
               py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> arr) {
                vector_assign_from_numpy(self.an, arr);
            },
            R"pbdoc(
                Returns the 'an' scattering coefficients.

                .. math::
                    a_n = \frac{
                        \Psi_n (y)
                        \left[ \dot{\Psi}_n (\tilde{m}_{2} y) - A_n (x) \dot{\xi}_n (\tilde{m}_{2} y) \right]
                        -
                        \tilde{m}_{2} \dot{\Psi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - A_n (x) \xi_n (\tilde{m}_{2} y) \right]
                    }{
                        \xi_n (y)
                        \left[ \dot{\Psi}_n (\tilde{m}_{2} y) - A_n (x) \dot{\xi}_n (\tilde{m}_{2} y) \right]
                        -
                        \tilde{m}_{2} \dot{\xi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - A_n (x) \xi_n (\tilde{m}_{2} y) \right]
                    }

                with:

                .. math::
                    A_n = \frac{
                        \tilde{m}_{2}  \Psi_n       (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\Psi}_n (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \xi_n        (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\xi}_n  (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }
                    \\
                    B_n = \frac{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\psi}_n (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \psi_n  (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\xi}_n  (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \xi_n   (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }

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
                Returns the 'bn' scattering coefficients.

                .. math::
                    b_n = \frac{
                        \tilde{m}_{2}
                        \Psi_n (y)
                        \left[ \dot{\Psi}_n  (\tilde{m}_{2} y) - B_n (x) \dot{\xi}_n  (\tilde{m}_{2} y) \right]
                        -
                        \dot{\Psi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - B_n (x) \xi_n  (\tilde{m}_{2} y) \right]
                    }{
                        \tilde{m}_{2}
                        \xi_n (y)
                        \left[ \dot{\Psi}_n  (\tilde{m}_{2} y) - A_n (x) \dot{\xi}_n  (\tilde{m}_{2} y) \right] -
                        \dot{\xi}_n (y)
                        \left[ \Psi_n (\tilde{m}_{2} y) - B_n (x) \xi_n (\tilde{m}_{2} y) \right]
                    }

                with:

                .. math::
                    A_n = \frac{
                        \tilde{m}_{2}  \Psi_n       (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\Psi}_n (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \xi_n        (\tilde{m}_{2} x) \dot{\Psi}_n (\tilde{m}_{1} x) -
                        \tilde{m}_{1}  \dot{\xi}_n  (\tilde{m}_{2} x) \Psi_n       (\tilde{m}_{1} x)
                    }
                    \\
                    B_n = \frac{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\psi}_n (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \psi_n  (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }{
                        \tilde{m}_{2}  \psi_n  (\tilde{m}_{1} x) \dot{\xi}_n  (\tilde{m}_{2} x) -
                        \tilde{m}_{1}  \xi_n   (\tilde{m}_{2} x) \dot{\psi}_n (\tilde{m}_{1} x)
                    }

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
                Returns the 'cn' scattering coefficients.

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
                Returns the 'dn' scattering coefficients.

                Returns
                -------
                list
                    A list of 'dn' scattering coefficients used in the spherical wave expansion.
            )pbdoc"
        )
        ;
}
