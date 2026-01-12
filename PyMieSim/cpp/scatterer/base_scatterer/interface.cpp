#include <pybind11/pybind11.h>
#include "base_scatterer.h"
#include "../../utils/numpy_interface.h"
#include "pint/pint.h"


void register_base_scatterer(pybind11::module_& module) {

    pybind11::class_<BaseScatterer>(module, "BASESCATTERER")
        .def_readonly("_cpp_cross_section",
            &BaseScatterer::cross_section,
            R"pbdoc(
                Cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The cross-sectional area (cross_section) used in scattering calculations.
            )pbdoc"
        )
        .def("_cpp_compute_cn_dn",
            &BaseScatterer::compute_cn_dn,
            pybind11::arg("max_order") = 0,
            R"pbdoc(
                Computes the coefficients cn and dn for the scatterer.

                Parameters
                ----------
                max_order : int, optional
                    The maximum order of the coefficients to compute (default is 0, which means it will be computed based on the size parameter).
            )pbdoc"
        )
        .def("_cpp_compute_an_bn",
            &BaseScatterer::compute_an_bn,
            pybind11::arg("max_order") = 0,
            R"pbdoc(
                Computes the coefficients an and bn for the scatterer.
                Parameters
                ----------
                max_order : int, optional
                    The maximum order of the coefficients to compute (default is 0, which means it will be computed based on the size parameter).
            )pbdoc"
        )
        .def("_cpp_get_s1s2",
            [](BaseScatterer& self, const std::vector<double> &phi){
                auto [S1, S2] = self.compute_s1s2(phi);
                return std::make_tuple(
                    vector_move_from_numpy(S1, {S1.size()}),
                    vector_move_from_numpy(S2, {S2.size()})
                );
            },
            pybind11::arg("phi"),
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
            )pbdoc"
        )
        .def("_cpp_get_farfields",
            [](BaseScatterer& self, const std::vector<double>& phi, const std::vector<double>& theta, const double& distance) {
                auto [theta_field, phi_field] = self.compute_unstructured_farfields(phi, theta, distance);

                return std::make_tuple(
                    vector_move_from_numpy(std::move(phi_field), {phi_field.size()}),
                    vector_move_from_numpy(std::move(theta_field), {theta_field.size()})
                );
            },
            pybind11::arg("phi"),
            pybind11::arg("theta"),
            pybind11::arg("distance"),
            pybind11::return_value_policy::move,
            R"pbdoc(
                Returns the unstructured electromagnetic fields around the scatterer.

                Parameters
                ----------
                phi : float
                    The azimuthal angle (in radians).
                theta : float
                    The polar angle (in radians).
                distance : float
                    The radial distance from the scatterer.

                Returns
                -------
                numpy.ndarray
                    A NumPy array containing the unstructured electromagnetic fields (E and H) around the scatterer.
            )pbdoc"
        )
        .def("_cpp_get_full_farfields",
            [](BaseScatterer& self, const size_t sampling, const double& distance) {

                auto [phi_field_mesh, theta_field_mesh, theta_mesh, phi_mesh] = self.compute_full_structured_farfields(sampling, distance);
                return std::make_tuple(
                    vector_move_from_numpy(std::move(phi_field_mesh), {sampling, sampling}).attr("transpose")(),
                    vector_move_from_numpy(std::move(theta_field_mesh), {sampling, sampling}).attr("transpose")(),
                    vector_move_from_numpy(std::move(phi_mesh), {sampling}),
                    vector_move_from_numpy(std::move(theta_mesh), {sampling})
                );
            },
            pybind11::arg("sampling"),
            pybind11::arg("distance"),
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
            )pbdoc"
        )
        .def("_cpp_get_coefficient",
            &BaseScatterer::get_coefficient,
            pybind11::arg("type"),
            pybind11::arg("order"),
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
            )pbdoc"
        )
        .def_property_readonly("Qsca",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_Qsca()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Scattering efficiency of the scatterer.

                .. math::
                    Q_{\text{sca}} &= \frac{2}{x^2} \sum_{n=1}^{n_{\text{max}}} (2n+1)(|a_n|^2 + |b_n|^2)

                Returns
                -------
                float
                    The scattering efficiency (Qsca) value, which characterizes the effectiveness of the scatterer in scattering incident light.
            )pbdoc"
        )
        .def_property_readonly("Qext",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_Qext()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Extinction efficiency of the scatterer.

                .. math::
                    Q_{\text{ext}} &= \frac{2}{x^2} \sum_{n=1}^{n_{\text{max}}} (2n+1) \operatorname{Re}(a_n + b_n)

                Returns
                -------
                float
                    The extinction efficiency (Qext) value, representing the combined effect of scattering and absorption.
            )pbdoc"
        )
        .def_property_readonly("Qabs",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_Qabs()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Absorption efficiency of the scatterer.

                .. math::
                    Q_{\text{abs}} &= Q_{\text{ext}} - Q_{\text{sca}}

                Returns
                -------
                float
                    The absorption efficiency (Qabs) value, which characterizes the absorption of incident light by the scatterer.
            )pbdoc"
        )
        .def_property_readonly("Qback",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_Qback()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Backscattering efficiency of the scatterer.

                .. math::
                    Q_{\text{back}} &= \frac{1}{x^2} \left| \sum_{n=1}^{n_{\text{max}}} (2n+1)(-1)^n (a_n - b_n) \right|^2

                Returns
                -------
                float
                    The backscattering efficiency (Qback) value, which characterizes the amount of light scattered in the backward direction.
            )pbdoc"
        )
        .def_property_readonly("Qforward",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_Qforward()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Forward-scattering efficiency of the scatterer.

                Returns
                -------
                float
                    The forward-scattering efficiency (Qforward) value, which characterizes the amount of light scattered in the forward direction.
            )pbdoc"
        )
        .def_property_readonly("Qratio",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_Qratio()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                The ratio of forward to backward scattering efficiency (Qforward/Qback) of the scatterer.

                .. math::
                    Q_{\text{ratio}} &= \frac{Q_{\text{back}}}{Q_{\text{sca}}}

                Returns
                -------
                float
                    The ratio of forward to backward scattering efficiency (Qratio), giving insight into the asymmetry of scattering.
            )pbdoc"
        )
        .def_property_readonly("Qpr",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_Qpr()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Radiation pressure efficiency of the scatterer.

                .. math::
                    Q_{\text{pr}} &= Q_{\text{ext}} - g \cdot Q_{\text{sca}}

                Returns
                -------
                float
                    The radiation pressure efficiency (Qpr) value, associated with the pressure exerted by the scattered light on the scatterer.
            )pbdoc"
        )
        .def_property_readonly("Csca",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.get_Csca()) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Scattering cross-section of the scatterer.

                .. math::
                    C_{\text{sca}} &= Q_{\text{sca}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The scattering cross-section (Csca) value, which represents the total scattering area of the scatterer.
            )pbdoc"
        )
        .def_property_readonly("Cext",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.get_Cext()) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Extinction cross-section of the scatterer.

                .. math::
                    C_{\text{ext}} &= Q_{\text{ext}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The extinction cross-section (Cext) value, representing the total extinction, including scattering and absorption.
            )pbdoc"
        )
        .def_property_readonly("Cabs",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.get_Cabs()) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Absorption cross-section of the scatterer.

                .. math::
                    C_{\text{abs}} &= Q_{\text{abs}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The absorption cross-section (Cabs) value, representing the area over which the scatterer absorbs light.
            )pbdoc"
        )
        .def_property_readonly("Cback",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.get_Cback()) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Backscattering cross-section of the scatterer.

                .. math::
                    C_{\text{back}} &= Q_{\text{back}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The backscattering cross-section (Cback) value, representing the total scattering area for backward-scattered light.
            )pbdoc"
        )
        .def_property_readonly("Cforward",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.get_Cforward()) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Forward-scattering cross-section of the scatterer.

                .. math::
                    C_{\text{forward}} &= Q_{\text{forward}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The forward-scattering cross-section (Cforward) value, representing the total scattering area for forward-scattered light.
            )pbdoc"
        )
        .def_property_readonly("Cratio",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.get_Cratio()) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                The ratio of forward to backward scattering cross-section (Cforward/Cback) of the scatterer.

                .. math::
                    C_{\text{ratio}} &= \frac{C_{\text{forward}}}{C_{\text{back}}}

                where Cforward is the forward scattering cross-section and Cback is the backward scattering cross-section.

                Returns
                -------
                float
                    The ratio of forward to backward scattering cross-section (Cratio), which gives an indication of scattering asymmetry.
            )pbdoc"
        )
        .def_property_readonly("Cpr",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.get_Cpr()) * ureg.attr("meter**2")).attr("to_compact")();;
            },
            R"pbdoc(
                Radiation pressure cross-section of the scatterer.

                .. math::
                    C_{\text{pr}} &= Q_{\text{pr}} \cdot \text{area}

                where area is the physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The radiation pressure cross-section (Cpr) value, associated with the cross-sectional area for radiation pressure exerted on the scatterer.
            )pbdoc"
        )
        .def_property_readonly("g",
            [](BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.get_g()) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Asymmetry parameter of the scatterer.

                .. math::
                    g &= \frac{4}{Q_{\text{sca}} x^2} \left[ \sum_{n=1}^{n_{\text{max}}} \frac{n(n+2)}{n+1} \operatorname{Re}(a_n a_{n+1}^* + b_n b_{n+1}^*) + \sum_{n=1}^{n_{\text{max}}} \frac{2n+1}{n(n+1)} \operatorname{Re}(a_n b_n^*) \right]


                Returns
                -------
                float
                    The asymmetry parameter (g) value, where g = 0 represents isotropic scattering.
            )pbdoc"
        )
        .def_property_readonly("cross_section",
            [](const BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return (py::float_(self.cross_section) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Physical cross-sectional area of the scatterer.

                Returns
                -------
                float
                    The physical cross-sectional area (area) used in scattering calculations.
            )pbdoc"
        )
        .def_property_readonly("size_parameter",
            [](const BaseScatterer& self) {
                py::object ureg = get_shared_ureg();
                return py::float_(self.size_parameter) * ureg.attr("dimensionless");
            },
            R"pbdoc(
                Size parameter of the scatterer.

                The size parameter is defined as the ratio of the scatterer's physical cross-section to the wavelength of incident light, scaled by the refractive index of the medium.

                Returns
                -------
                float
                    The size parameter (size_parameter), typically the ratio of the scatterer's diameter to the wavelength of incident light.
            )pbdoc"
        )
        .def("_cpp_compute_nearfields",  // &BaseScatterer::compute_nearfields_py,
            [](BaseScatterer& self, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::string& field_type) {
                std::vector<complex128> field = self.compute_nearfields(x, y, z, field_type);

                return vector_move_from_numpy(field, {field.size()});
            },
            pybind11::arg("x"),
            pybind11::arg("y"),
            pybind11::arg("z"),
            pybind11::arg("field_type"),
            R"pbdoc(
                Computes near-field electromagnetic fields using internal coefficients cn and dn.

                This method calculates the electromagnetic fields inside and near the scatterer
                using the multipole expansion with vector spherical harmonics. The internal fields
                (r < radius) are computed using cn and dn coefficients, while external fields
                (r > radius) use an and bn coefficients.

                Parameters
                ----------
                x : numpy.ndarray
                    Array of x coordinates of observation points.
                y : numpy.ndarray
                    Array of y coordinates of observation points.
                z : numpy.ndarray
                    Array of z coordinates of observation points.
                field_type : str
                    Field component type: "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "|E|", "|H|"

                Returns
                -------
                numpy.ndarray
                    Complex array of field values at specified points.

                Raises
                ------
                RuntimeError
                    If cn/dn coefficients are not available for the scatterer type.
                ValueError
                    If field_type is not recognized.

                Notes
                -----
                This method requires that cn and dn coefficients have been computed.
                Currently supports spherical scatterers only. Cylinder support requires
                implementation of cn/dn coefficients for infinite cylinders.
            )pbdoc"
        )
        .def("_cpp_compute_nearfields_structured",
            [](BaseScatterer& self, const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z, const std::string& field_type)
            {
                auto [field_values, field_x, field_y, field_z, x_coords, y_coords, z_coords] = self.compute_nearfields_structured(x, y, z, field_type);

                return std::make_tuple(
                    vector_move_from_numpy(field_values, {x.size(), y.size(), z.size()}),
                    vector_move_from_numpy(field_x, {x.size(), y.size(), z.size()}),
                    vector_move_from_numpy(field_y, {x.size(), y.size(), z.size()}),
                    vector_move_from_numpy(field_z, {x.size(), y.size(), z.size()}),
                    vector_move_from_numpy(x_coords, {x.size() * y.size() * z.size()}),
                    vector_move_from_numpy(y_coords, {x.size() * y.size() * z.size()}),
                    vector_move_from_numpy(z_coords, {x.size() * y.size() * z.size()})
                );
            },
            pybind11::arg("x"),
            pybind11::arg("y"),
            pybind11::arg("z"),
            pybind11::arg("field_type"),
            R"pbdoc(
                Computes near-field electromagnetic fields over a structured 3D grid.

                This method efficiently computes electromagnetic fields over a regular 3D grid
                using the multipole expansion with vector spherical harmonics. More efficient
                than point-by-point computation for regular grids.

                Parameters
                ----------
                x : numpy.ndarray
                    Array of x coordinates defining the grid.
                y : numpy.ndarray
                    Array of y coordinates defining the grid.
                z : numpy.ndarray
                    Array of z coordinates defining the grid.
                field_type : str
                    Field component type: "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "|E|", "|H|"

                Returns
                -------
                tuple
                    Tuple containing:
                    - Requested field component values (nx . ny . nz) as numpy.ndarray
                    - Ex field component (nx . ny . nz) as numpy.ndarray
                    - Ey field component (nx . ny . nz) as numpy.ndarray
                    - Ez field component (nx . ny . nz) as numpy.ndarray
                    - x coordinates (flattened) as numpy.ndarray
                    - y coordinates (flattened) as numpy.ndarray
                    - z coordinates (flattened) as numpy.ndarray

                Raises
                ------
                RuntimeError
                    If cn/dn coefficients are not available.
                ValueError
                    If field_type is not recognized.

                Notes
                -----
                This method provides all field components for comprehensive analysis and
                visualization of near-field effects around the scatterer.
            )pbdoc"
        )
        ;
}
