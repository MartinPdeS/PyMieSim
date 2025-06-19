#include <pybind11/pybind11.h>
#include "base_scatterer.h"


void register_base_scatterer(pybind11::module_& module) {

    pybind11::class_<BaseScatterer>(module, "BASESCATTERER")
        .def_readonly("_cpp_cross_section", &BaseScatterer::cross_section,
            R"pbdoc(
                Cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The cross-sectional area (cross_section) used in scattering calculations.
            )pbdoc")
        .def("_cpp_get_s1s2", &BaseScatterer::get_s1s2_py, pybind11::arg("phi"),
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
            pybind11::arg("phi"), pybind11::arg("theta"), pybind11::arg("r"), pybind11::return_value_policy::move,
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
            pybind11::arg("sampling"), pybind11::arg("distance"),
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
            pybind11::arg("type"), pybind11::arg("order"),
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
}

