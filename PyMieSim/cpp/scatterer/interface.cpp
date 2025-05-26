#include <pybind11/pybind11.h>
#include <pybind11/complex.h> // For complex number support
#include <pybind11/numpy.h>
#include "scatterer/sphere/sphere.h"
#include "scatterer/coreshell/coreshell.h"
#include "scatterer/cylinder/cylinder.h"



namespace py = pybind11;

PYBIND11_MODULE(interface_scatterer, module) {
    module.doc() = R"pbdoc(
        Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.

        This module provides C++ bindings for the PyMieSim Python package, which implements the Lorenz-Mie Theory (LMT) for light scattering by spherical particles and other scatterers.
    )pbdoc";

    py::class_<BaseScatterer>(module, "BASESCATTERER")
        .def_readonly("_cpp_cross_section", &BaseScatterer::cross_section,
            R"pbdoc(
                Cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The cross-sectional area (cross_section) used in scattering calculations.
            )pbdoc")
        .def("_cpp_get_s1s2", &BaseScatterer::get_s1s2_py, py::arg("phi"),
            R"pbdoc(
                Calculates and returns the S1 and S2 scattering parameters.

                Parameters
                ----------
                phi : float
                    The angle (in radians) at which the scattering is calculated.

                Returns
                -------
                tuple
                    A tuple containing the S1 and S2 scattering parameters, which represent the scattering amplitudes for the incident and scattered waves.
            )pbdoc")
        .def("_cpp_get_fields", &BaseScatterer::get_unstructured_fields_py,
            py::arg("phi"), py::arg("theta"), py::arg("r"), py::return_value_policy::move,
            R"pbdoc(
                Returns the unstructured electromagnetic fields around the scatterer.

                Parameters
                ----------
                phi : float
                    The azimuthal angle (in radians).
                theta : float
                    The polar angle (in radians).
                r : float
                    The radial distance from the scatterer.

                Returns
                -------
                numpy.ndarray
                    A NumPy array containing the unstructured electromagnetic fields (E and H) around the scatterer.
            )pbdoc")
        .def("_cpp_get_full_fields", &BaseScatterer::get_full_structured_fields_py,
            py::arg("sampling"), py::arg("distance"),
            R"pbdoc(
                Returns the full structured electromagnetic fields around the scatterer.

                Parameters
                ----------
                sampling : float
                    The sampling rate used in the field calculation.
                distance : float
                    The radial distance from the scatterer.

                Returns
                -------
                numpy.ndarray
                    A NumPy array with the full structured electromagnetic fields around the scatterer.
            )pbdoc")
        .def("_cpp_get_coefficient", &BaseScatterer::get_coefficient_py,
            py::arg("type"), py::arg("order"),
            R"pbdoc(
                Returns a specific scattering coefficient for a given type and order.

                Parameters
                ----------
                type : str
                    The type of coefficient ('an', 'bn', 'cn', 'dn').
                order : int
                    The order of the scattering coefficient.

                Returns
                -------
                float
                    The value of the specified scattering coefficient.
            )pbdoc")
        .def_property_readonly("_cpp_Qsca", &BaseScatterer::get_Qsca,
            R"pbdoc(
                Scattering efficiency of the scatterer.

                .. math::
                    Q_{\text{sca}} &= \frac{2}{x^2} \sum_{n=1}^{n_{\text{max}}} (2n+1)(|a_n|^2 + |b_n|^2)

                Returns
                -------
                float
                    The scattering efficiency (Qsca) value, which characterizes the effectiveness of the scatterer in scattering incident light.
            )pbdoc")
        .def_property_readonly("_cpp_Qext", &BaseScatterer::get_Qext,
            R"pbdoc(
                Extinction efficiency of the scatterer.

                .. math::
                    Q_{\text{ext}} &= \frac{2}{x^2} \sum_{n=1}^{n_{\text{max}}} (2n+1) \operatorname{Re}(a_n + b_n)

                Returns
                -------
                float
                    The extinction efficiency (Qext) value, representing the combined effect of scattering and absorption.
            )pbdoc")
        .def_property_readonly("_cpp_Qabs", &BaseScatterer::get_Qabs,
            R"pbdoc(
                Absorption efficiency of the scatterer.

                .. math::
                    Q_{\text{abs}} &= Q_{\text{ext}} - Q_{\text{sca}}

                Returns
                -------
                float
                    The absorption efficiency (Qabs) value, which characterizes the absorption of incident light by the scatterer.
            )pbdoc")
        .def_property_readonly("_cpp_Qback", &BaseScatterer::get_Qback,
            R"pbdoc(
                Backscattering efficiency of the scatterer.

                .. math::
                    Q_{\text{back}} &= \frac{1}{x^2} \left| \sum_{n=1}^{n_{\text{max}}} (2n+1)(-1)^n (a_n - b_n) \right|^2

                Returns
                -------
                float
                    The backscattering efficiency (Qback) value, which characterizes the amount of light scattered in the backward direction.
            )pbdoc")
        .def_property_readonly("_cpp_Qforward", &BaseScatterer::get_Qforward,
            R"pbdoc(
                Forward-scattering efficiency of the scatterer.

                Returns
                -------
                float
                    The forward-scattering efficiency (Qforward) value, which characterizes the amount of light scattered in the forward direction.
            )pbdoc")
        .def_property_readonly("_cpp_Qratio", &BaseScatterer::get_Qratio,
            R"pbdoc(
                The ratio of forward to backward scattering efficiency (Qforward/Qback) of the scatterer.

                .. math::
                    Q_{\text{ratio}} &= \frac{Q_{\text{back}}}{Q_{\text{sca}}}

                Returns
                -------
                float
                    The ratio of forward to backward scattering efficiency (Qratio), giving insight into the asymmetry of scattering.
            )pbdoc")
        .def_property_readonly("_cpp_Qpr", &BaseScatterer::get_Qpr,
            R"pbdoc(
                Radiation pressure efficiency of the scatterer.

                .. math::
                    Q_{\text{pr}} &= Q_{\text{ext}} - g \cdot Q_{\text{sca}}

                Returns
                -------
                float
                    The radiation pressure efficiency (Qpr) value, associated with the pressure exerted by the scattered light on the scatterer.
            )pbdoc")
        .def_property_readonly("_cpp_Csca", &BaseScatterer::get_Csca,
            R"pbdoc(
                Scattering cross-section of the scatterer.

                .. math::
                    C_{\text{sca}} &= Q_{\text{sca}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The scattering cross-section (Csca) value, which represents the total scattering area of the scatterer.
            )pbdoc")
        .def_property_readonly("_cpp_Cext", &BaseScatterer::get_Cext,
            R"pbdoc(
                Extinction cross-section of the scatterer.

                .. math::
                    C_{\text{ext}} &= Q_{\text{ext}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The extinction cross-section (Cext) value, representing the total extinction, including scattering and absorption.
            )pbdoc")
        .def_property_readonly("_cpp_Cabs", &BaseScatterer::get_Cabs,
            R"pbdoc(
                Absorption cross-section of the scatterer.

                .. math::
                    C_{\text{abs}} &= Q_{\text{abs}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The absorption cross-section (Cabs) value, representing the area over which the scatterer absorbs light.
            )pbdoc")
        .def_property_readonly("_cpp_Cback", &BaseScatterer::get_Cback,
            R"pbdoc(
                Backscattering cross-section of the scatterer.

                .. math::
                    C_{\text{back}} &= Q_{\text{back}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The backscattering cross-section (Cback) value, representing the total scattering area for backward-scattered light.
            )pbdoc")
        .def_property_readonly("_cpp_Cforward", &BaseScatterer::get_Cforward,
            R"pbdoc(
                Forward-scattering cross-section of the scatterer.

                .. math::
                    C_{\text{forward}} &= Q_{\text{forward}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The forward-scattering cross-section (Cforward) value, representing the total scattering area for forward-scattered light.
            )pbdoc")
        .def_property_readonly("_cpp_Cratio", &BaseScatterer::get_Cratio,
            R"pbdoc(
                The ratio of forward to backward scattering cross-section (Cforward/Cback) of the scatterer.

                .. math::
                    C_{\text{ratio}} &= \frac{C_{\text{forward}}}{C_{\text{back}}}

                where Cforward is the forward scattering cross-section and Cback is the backward scattering cross-section.

                Returns
                -------
                float
                    The ratio of forward to backward scattering cross-section (Cratio), which gives an indication of scattering asymmetry.
            )pbdoc")
        .def_property_readonly("_cpp_Cpr", &BaseScatterer::get_Cpr,
            R"pbdoc(
                Radiation pressure cross-section of the scatterer.

                .. math::
                    C_{\text{pr}} &= Q_{\text{pr}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The radiation pressure cross-section (Cpr) value, associated with the cross-sectional area for radiation pressure exerted on the scatterer.
            )pbdoc")
        .def_property_readonly("_cpp_g", &BaseScatterer::get_g,
            R"pbdoc(
                Asymmetry parameter of the scatterer.

                .. math::
                    g &= \frac{4}{Q_{\text{sca}} x^2} \left[ \sum_{n=1}^{n_{\text{max}}} \frac{n(n+2)}{n+1} \operatorname{Re}(a_n a_{n+1}^* + b_n b_{n+1}^*) + \sum_{n=1}^{n_{\text{max}}} \frac{2n+1}{n(n+1)} \operatorname{Re}(a_n b_n^*) \right]


                Returns
                -------
                float
                    The asymmetry parameter (g) value, where g = 0 represents isotropic scattering.
            )pbdoc")
        .def_readwrite("_cpp_cross_section", &BaseScatterer::cross_section,
            R"pbdoc(
                Physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The physical cross-sectional area (area) used in scattering calculations.
            )pbdoc")
        .def_readwrite("_cpp_size_parameter", &BaseScatterer::size_parameter,
            R"pbdoc(
                Size parameter of the scatterer.

                The size parameter is defined as the ratio of the scatterer's physical cross-section to the wavelength of incident light, scaled by the refractive index of the medium.

                Returns
                -------
                float
                    The size parameter (size_parameter), typically the ratio of the scatterer's diameter to the wavelength of incident light.
            )pbdoc")
        ;


    // Binding for Sphere class
    py::class_<Sphere, BaseScatterer>(module, "SPHERE")
        .def(
            py::init<const double, const complex128, const double, const BaseSource&, size_t>(),
            py::arg("diameter"),
            py::arg("refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
            py::arg("max_order") = 0,
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
        .def("an", &BaseScatterer::get_an_list_py,
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
            )pbdoc")
        .def("bn", &BaseScatterer::get_bn_list_py,
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
            )pbdoc")
        .def("cn", &BaseScatterer::get_cn_list_py,
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
            )pbdoc")
        .def("dn", &BaseScatterer::get_dn_list_py,
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

    // Binding for CoreShell class
    py::class_<CoreShell, BaseScatterer>(module, "CORESHELL")
        .def(py::init<double, double, std::complex<double>, std::complex<double>, double, BaseSource&>(),
            py::arg("core_diameter"),
            py::arg("shell_thickness"),
            py::arg("core_refractive_index"),
            py::arg("shell_refractive_index"),
            py::arg("medium_refractive_index"),
            py::arg("source"),
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
        .def("an", &BaseScatterer::get_an_list_py,
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
            )pbdoc")
        .def("bn", &BaseScatterer::get_bn_list_py,
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
            )pbdoc")
        .def("cn", &BaseScatterer::get_cn_list_py,
            R"pbdoc(
                Returns the 'cn' scattering coefficients.

                Returns
                -------
                list
                    A list of 'cn' scattering coefficients used in the spherical wave expansion.
            )pbdoc")
        .def("dn", &BaseScatterer::get_dn_list_py,
            R"pbdoc(
                Returns the 'dn' scattering coefficients.

                Returns
                -------
                list
                    A list of 'dn' scattering coefficients used in the spherical wave expansion.
            )pbdoc")
        ;

    // Binding for Cylinder class
    py::class_<Cylinder, BaseScatterer>(module, "CYLINDER")
        .def(
            py::init<double, complex128, double, BaseSource&>(),
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
        .def("a1n", &Cylinder::get_a1n_list_py,
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
        .def("b1n", &Cylinder::get_b1n_list_py,
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
        .def("a2n", &Cylinder::get_a2n_list_py,
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
        .def("b2n", &Cylinder::get_b2n_list_py,
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

