#include <pybind11/pybind11.h>

#include <utils/numpy_interface.h>
#include <utils/constants.h>
#include "./setup.h"
#include <pint/pint.h>

namespace py = pybind11;

#define RETURN_PROPERTY_UNITS(name, units) \
    if (property == #name) return py::float_(self.scatterer->get_##name()) * ureg.attr(#units);

#define RETURN_PROPERTY(name) \
    if (property == #name) return py::float_(self.scatterer->get_##name()) * ureg.attr("dimensionless");


PYBIND11_MODULE(setup, module)
{
    py::object ureg = get_shared_ureg();
    py::module_::import("PyMieSim.single.detector");
    py::module_::import("PyMieSim.single.source");
    py::module_::import("PyMieSim.single.scatterer");
    py::module_::import("PyMieSim.mesh");

    py::class_<Setup, std::shared_ptr<Setup>>(module,
        "Setup",
        py::module_local(),
            R"pbdoc(
                Setup for single scatterer Lorenz-Mie simulations.

                This class encapsulates the scatterer, source, and detector
                configurations for a single scattering simulation.
            )pbdoc"
        )
        .def_readonly(
            "source",
            &Setup::source,
            R"pbdoc(
                The source configuration for the simulation.
            )pbdoc"
        )
        .def_readonly(
            "scatterer",
            &Setup::scatterer,
            R"pbdoc(
                The scatterer configuration for the simulation.
            )pbdoc"
        )
        .def_readonly(
            "detector",
            &Setup::detector,
            R"pbdoc(
                The detector configuration for the simulation.
            )pbdoc"
        )
        .def(
            py::init(
                [](
                    std::shared_ptr<BaseScatterer> scatterer,
                    std::shared_ptr<BaseSource> source,
                    std::shared_ptr<BaseDetector> detector = nullptr,
                    bool debug_mode = false
                ) {
                    return std::make_shared<Setup>(scatterer, source, detector, debug_mode);
                }
            ),
            py::arg("scatterer"),
            py::arg("source"),
            py::arg("detector") = nullptr,
            py::arg("debug_mode") = false,
            R"pbdoc(
                Initializes a Setup instance for single scatterer simulations.
                Parameters
                ----------
                scatterer : BaseScatterer
                    The scatterer configuration.
                source : BaseSource
                    The source configuration.
                detector : Photodiode
                    The detector configuration.
                debug_mode : bool, optional
                    If True, enables debug mode with additional logging. Default is False.
            )pbdoc"
        )
        .def("get",
            [ureg](const Setup& self, const std::string& property) {
                RETURN_PROPERTY(Qsca)
                RETURN_PROPERTY_UNITS(Csca, meter**2)

                RETURN_PROPERTY(Qext)
                RETURN_PROPERTY_UNITS(Cext, meter**2)

                RETURN_PROPERTY(Qabs)
                RETURN_PROPERTY_UNITS(Cabs, meter**2)

                RETURN_PROPERTY(Qback)
                RETURN_PROPERTY_UNITS(Cback, meter**2)

                RETURN_PROPERTY(Qforward)
                RETURN_PROPERTY_UNITS(Cforward, meter**2)

                RETURN_PROPERTY(Qratio)
                RETURN_PROPERTY_UNITS(Cratio, meter**2)

                RETURN_PROPERTY(Qpr)
                RETURN_PROPERTY_UNITS(Cpr, meter**2)

                RETURN_PROPERTY(g)
                if (property == "g_with_farfields") return py::float_(self.scatterer->get_g_with_farfields(self.source, 1000)) * ureg.attr("dimensionless");

                RETURN_PROPERTY_UNITS(cross_section, meter**2)
                RETURN_PROPERTY(size_parameter)

                if (property == "coupling") {
                    if (!self.detector) {
                        throw std::logic_error("Detector is not defined in the setup. Coupling cannot be computed.");
                    }
                    return py::float_(self.detector->get_coupling(self.scatterer, self.source)) * ureg.attr("watt");

                }

                throw std::invalid_argument("Unknown data name: " + property);
            },
            py::arg("data"),
            R"pbdoc(
                Retrieves specific data from the simulation results based on the provided data name.

                Parameters
                ----------
                data_name : str
                    The name of the data to retrieve. Supported values include 'coupling', 'Qsca', 'Qext', 'Qback', 'Qforward', and 'g'.

                    Returns
                -------
                float
                    The requested data value corresponding to the provided data name.

                Raises
                ------
                ValueError
                    If an unknown data name is provided.
            )pbdoc"
        )

        //------------------------ S1S2 --------------------------

        .def(
            "get_s1s2",
            [ureg](Setup& self, const py::object& angles)
            {
                std::vector<double> angle_values =
                    angles.attr("to")("radian")
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                std::pair<std::vector<complex128>, std::vector<complex128>> s1s2 =
                    self.get_s1s2(angle_values);

                py::array_t<complex128> s1_array = _vector_to_numpy(s1s2.first);
                py::array_t<complex128> s2_array = _vector_to_numpy(s1s2.second);

                return py::make_tuple(
                    s1_array * ureg.attr("meter"),
                    s2_array * ureg.attr("meter")
                );
            },
            py::arg("angles"),
            R"pbdoc(
                Compute the complex scattering amplitude functions :math:`S_1(\phi)` and :math:`S_2(\phi)`.

                These quantities are the standard amplitude scattering functions used in Lorenz Mie theory.

                They relate the incident field to the scattered far field and define the angular dependence of the scattered response for the two orthogonal polarization channels.

                In this implementation, :math:`S_1` and :math:`S_2` are returned with units of length.

                The differential scattering cross section is obtained from these amplitudes as

                .. math::

                    \frac{d\sigma}{d\Omega} = |S_1(\phi)|^2

                or

                .. math::

                    \frac{d\sigma}{d\Omega} = |S_2(\phi)|^2

                depending on the incident polarization basis, and for unpolarized illumination one commonly writes

                .. math::

                    \frac{d\sigma}{d\Omega} = \frac{|S_1(\phi)|^2 + |S_2(\phi)|^2}{2}.

                Parameters
                ----------
                angles : pint.Quantity
                    One dimensional array of scattering angles. Values must be convertible to radians.

                Returns
                -------
                tuple[pint.Quantity, pint.Quantity]
                    Two one dimensional arrays ``(S1, S2)`` with complex dtype and units of meter.

                Notes
                -----
                The returned arrays have the same length and ordering as the input ``angles``.
                %
                The angle convention used here is the internal convention of the scattering solver.
                %
                In the underlying C++ implementation, the input angles are shifted by :math:`\pi/2` before calling the scatterer level routine, so this method should be treated as the public high level angular interface for users.
            )pbdoc"
        )



        //------------------------ STOKES --------------------------
        .def(
            "get_stokes",
            [ureg](
                Setup& self,
                const size_t sampling,
                const py::object& distance
            ) {
                std::tuple<
                    std::vector<double>,
                    std::vector<double>,
                    std::vector<double>,
                    std::vector<double>,
                    FullSteradian
                > stokes_parameters = // [i, Q, U, V, FullSteradian]
                     self.get_structured_stokes(
                        sampling,
                        distance.attr("to")("meter").attr("magnitude").cast<double>()
                );

                std::vector<size_t> shape = {sampling, sampling};
                py::array_t<double> i_array = _vector_to_numpy(std::get<0>(stokes_parameters), shape);
                py::array_t<double> Q_array = _vector_to_numpy(std::get<1>(stokes_parameters), shape);
                py::array_t<double> U_array = _vector_to_numpy(std::get<2>(stokes_parameters), shape);
                py::array_t<double> V_array = _vector_to_numpy(std::get<3>(stokes_parameters), shape);

                return py::make_tuple(
                    i_array * ureg.attr("watt/m**2"),
                    Q_array * ureg.attr("watt/m**2"),
                    U_array * ureg.attr("watt/m**2"),
                    V_array * ureg.attr("watt/m**2"),
                    std::get<4>(stokes_parameters)
                );
            },
            py::arg("sampling"),
            py::arg("distance"),
             R"pbdoc(
                Returns the unstructured Stokes parameters around the scatterer.

                Parameters
                ----------
                sampling : int
                    The number of points in the angular sampling for both phi and theta directions.
                distance : float
                    The radial distance from the scatterer.

                Returns
                -------
                tuple
                    A NumPy array with the unstructured Stokes parameters (I, Q, U, V) around the scatterer.

            )pbdoc"
        )
        .def(
            "get_stokes",
            [ureg](
                Setup& self,
                const py::object& phi,
                const py::object& theta,
                const py::object& distance
            ) {
                std::tuple<
                    std::vector<double>,
                    std::vector<double>,
                    std::vector<double>,
                    std::vector<double>
                > stokes_parameters = self.get_unstructured_stokes(
                    phi.attr("to")("radian").attr("magnitude").cast<std::vector<double>>(),
                    theta.attr("to")("radian").attr("magnitude").cast<std::vector<double>>(),
                    distance.attr("to")("meter").attr("magnitude").cast<double>()
                );

                py::array_t<double> i_array = _vector_to_numpy(std::get<0>(stokes_parameters));
                py::array_t<double> Q_array = _vector_to_numpy(std::get<1>(stokes_parameters));
                py::array_t<double> U_array = _vector_to_numpy(std::get<2>(stokes_parameters));
                py::array_t<double> V_array = _vector_to_numpy(std::get<3>(stokes_parameters));

                return py::make_tuple(
                    i_array * ureg.attr("watt/m**2"),
                    Q_array * ureg.attr("watt/m**2"),
                    U_array * ureg.attr("watt/m**2"),
                    V_array * ureg.attr("watt/m**2")
                );
            },
            py::arg("phi"),
            py::arg("theta"),
            py::arg("distance"),
            R"pbdoc(
                Returns the Stokes parameters at specified angles and distance from the scatterer.

                Parameters
                ----------
                phi : pint.Quantity
                    One dimensional array of azimuthal angles. Values must be convertible to radians.
                theta : pint.Quantity
                    One dimensional array of polar angles. Values must be convertible to radians.
                distance : pint.Quantity
                    Observation distance from the scatterer center. Values must be convertible to meter.

                Returns
                -------
                tuple
                    A tuple containing four one-dimensional NumPy arrays representing the Stokes parameters (I, Q, U, V) at the specified angles and distance from the scatterer.
            )pbdoc"
        )



        //------------------------ FARFIELDS --------------------------
        .def(
            "get_farfields",
            [ureg](
                Setup& self,
                const size_t sampling,
                const py::object& distance
            ) {
                std::tuple<
                    std::vector<complex128>,
                    std::vector<complex128>,
                    FullSteradian
                > farfield_meshes = self.get_structured_farfields(
                    sampling,
                    distance.attr("to")("meter").attr("magnitude").cast<double>()
                );

                std::vector<size_t> shape = {sampling, sampling};

                py::array_t<complex128> phi_field_mesh = _vector_to_numpy(std::get<0>(farfield_meshes), shape);
                py::array_t<complex128> theta_field_mesh = _vector_to_numpy(std::get<1>(farfield_meshes), shape);

                return py::make_tuple(
                    phi_field_mesh * ureg.attr("volt/meter"),
                    theta_field_mesh * ureg.attr("volt/meter"),
                    std::get<2>(farfield_meshes)
                );
            },
            py::arg("sampling"),
            py::arg("distance"),
            py::return_value_policy::copy,
            R"pbdoc(
                Computes the full structured far fields over a spherical mesh.

                Parameters
                ----------
                sampling : int
                    The number of points in the angular sampling for both phi and theta directions.
                distance : float
                    The radial distance from the scatterer at which to compute the far fields.

                Returns
                -------
                tuple
                    A tuple containing:
                    - phi_field_mesh: A 2D NumPy array of complex values representing the phi component of the far field on the mesh.
                    - theta_field_mesh: A 2D NumPy array of complex values representing the theta component of the far field on the mesh.
                    - phi_mesh: A 2D NumPy array of floats representing the phi angles of the mesh.
                    - theta_mesh: A 2D NumPy array of floats representing the theta angles of the mesh.
            )pbdoc"
        )
        .def(
            "get_farfields",
            [ureg](
                Setup& self,
                const py::object& phi,
                const py::object& theta,
                const py::object& distance
            ) {
                std::vector<double> phi_values =
                    phi.attr("to")("radian")
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                std::vector<double> theta_values =
                    theta.attr("to")("radian")
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                double distance_value =
                    distance.attr("to")("meter")
                    .attr("magnitude")
                    .cast<double>();

                std::pair<std::vector<complex128>, std::vector<complex128>> farfield_components =
                    self.get_unstructured_farfields(
                        phi_values,
                        theta_values,
                        distance_value
                    );

                py::array_t<complex128> e_phi_array = _vector_to_numpy(farfield_components.first);
                py::array_t<complex128> e_theta_array = _vector_to_numpy(farfield_components.second);

                return py::make_tuple(
                    e_phi_array * ureg.attr("volt/meter"),
                    e_theta_array * ureg.attr("volt/meter")
                );
            },
            py::arg("phi"),
            py::arg("theta"),
            py::arg("distance"),
            R"pbdoc(
                Compute the complex far field electric field components :math:`E_\phi` and :math:`E_\theta`.

                This method evaluates the scattered electric field in spherical coordinates at the requested angular positions and radial distance.
                %
                The returned quantities are the transverse far field components in the :math:`\phi` and :math:`\theta` directions.

                In the far field, the scattered electric field can be written in the generic form

                .. math::

                    \mathbf{E}_{\mathrm{sca}}(r, \theta, \phi)
                    \propto
                    \frac{e^{ikr}}{r}
                    \left[
                        E_\theta(\theta, \phi)\,\hat{\boldsymbol{\theta}}
                        +
                        E_\phi(\theta, \phi)\,\hat{\boldsymbol{\phi}}
                    \right],

                where :math:`r` is the observation distance and :math:`k` is the wavenumber in the surrounding medium.

                Parameters
                ----------
                phi : pint.Quantity
                    One dimensional array of azimuthal angles.
                    %
                    Values must be convertible to radians.

                theta : pint.Quantity
                    One dimensional array of polar angles.
                    %
                    Values must be convertible to radians.

                distance : pint.Quantity
                    Observation distance from the scatterer center.
                    %
                    Must be convertible to meter.

                Returns
                -------
                tuple[pint.Quantity, pint.Quantity]
                    Two one dimensional arrays ``(E_phi, E_theta)`` with complex dtype and units of volt per meter.

                Notes
                -----
                The arrays ``phi`` and ``theta`` are passed directly to the underlying far field solver.
                %
                The returned arrays have the same shape and ordering as the angular sampling defined by the input coordinates.
                %
                This method returns electric field components only. It does not return magnetic field components.
            )pbdoc"
        )




        //------------------------ SPF --------------------------
        .def(
            "get_spf",
            [ureg](
                Setup& self,
                const size_t sampling,
                const py::object& distance
            ) {
                std::pair<std::vector<double>, FullSteradian> spf_mesh = self.get_structured_spf(
                    sampling,
                    distance.attr("to")("meter").attr("magnitude").cast<double>()
                );

                std::vector<size_t> shape = {sampling, sampling};

                py::array_t<double> spf_field_mesh = _vector_to_numpy(std::get<0>(spf_mesh), shape);

                return py::make_tuple(
                    spf_field_mesh,
                    std::get<1>(spf_mesh)
                );
            },
            py::arg("sampling"),
            py::arg("distance") = py::float_(1.0) * ureg.attr("meter"),
            py::return_value_policy::copy,
            R"pbdoc(
                Computes the full structured scattering phase function (SPF) over a spherical mesh.

                Parameters
                ----------
                sampling : int
                    The number of points in the angular sampling for both phi and theta directions.
                distance : float, optional
                    The distance from the scatterer at which the scattering phase function is evaluated (default is 1.0 meter).

                Returns
                -------
                tuple
                    A tuple containing:
                    - spf_field_mesh: A 2D NumPy array of floats representing the scattering phase function values on the mesh.
                    - phi_mesh: A 2D NumPy array of floats representing the phi angles of the mesh.
                    - theta_mesh: A 2D NumPy array of floats representing the theta angles of the mesh.
            )pbdoc"
        )
        .def(
            "get_spf",
            [ureg](
                Setup& self,
                const py::object& phi,
                const py::object& theta,
                const py::object& distance
            ) {
                std::vector<double> phi_values =
                    phi.attr("to")("radian")
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                std::vector<double> theta_values =
                    theta.attr("to")("radian")
                    .attr("magnitude")
                    .cast<std::vector<double>>();

                double distance_value =
                    distance.attr("to")("meter")
                    .attr("magnitude")
                    .cast<double>();

                std::vector<double> spf_values = self.get_unstructured_spf(
                    phi_values,
                    theta_values,
                    distance_value
                );

                py::array_t<double> spf_array = _vector_to_numpy(spf_values);

                return spf_array;
            },
            py::arg("phi"),
            py::arg("theta"),
            py::arg("distance") = py::float_(1.0) * ureg.attr("meter"),
            R"pbdoc(
                Compute the unstructured scattering phase function (SPF) at specified angles and distance from the scatterer.

                Parameters
                ----------
                phi : pint.Quantity
                    One dimensional array of azimuthal angles. Values must be convertible to radians.

                theta : pint.Quantity
                    One dimensional array of polar angles. Values must be convertible to radians.

                distance : pint.Quantity
                    Observation distance from the scatterer center. Must be convertible to meter.

                Returns
                -------
                numpy.ndarray
                    A one-dimensional NumPy array of floats representing the scattering phase function values at the specified angles and distance from the scatterer.
            )pbdoc"
        )











        .def(
            "get_scattered_nearfields",
            [ureg](
                Setup& self,
                const py::object x,
                const py::object y,
                const py::object z,
                const std::string& field_type
            ) {
                std::vector<double> x_values = x.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();
                std::vector<double> y_values = y.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();
                std::vector<double> z_values = z.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                std::vector<complex128> neafields = self.get_scattered_nearfields(x_values, y_values, z_values, field_type);

                py::array_t<complex128> nearfield_array = _vector_to_numpy(neafields);

                return nearfield_array * ureg.attr("volt/meter");
            },
            py::arg("x"),
            py::arg("y"),
            py::arg("z"),
            py::arg("field_type"),
            R"pbdoc(
                Computes scattered near-field electromagnetic fields using an and bn coefficients.

                This method calculates the scattered electromagnetic fields outside the scatterer
                using the multipole expansion with vector spherical harmonics. It uses an and bn
                coefficients for the computation.

                Parameters
                ----------
                x : numpy.ndarray
                    Array of x coordinates of observation points.
                y : numpy.ndarray
                    Array of y coordinates of observation points.
                z : numpy.ndarray
                    Array of z coordinates of observation points.
                field_type : str
                    Field component type: "Ex", "Ey", "Ez", "\|E\|"
                source : BaseSource
                    The source of the electromagnetic field.

                Returns
                -------
                numpy.ndarray
                    Complex array of field values at specified points.

                Raises
                ------
                RuntimeError
                    If an/bn coefficients are not available for the scatterer type.
                ValueError
                    If field_type is not recognized.

                Notes
                -----
                This method requires that an and bn coefficients have been computed.
                Currently supports spherical scatterers only. Cylinder support requires
                implementation of an/bn coefficients for infinite cylinders.
            )pbdoc"
        )
        .def(
            "get_incident_nearfields",
            [ureg](
                Setup& self,
                const py::object x,
                const py::object y,
                const py::object z,
                const std::string& field_type
            ) {

                std::vector<complex128> field = self.get_incident_nearfields(
                    x.attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                    y.attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                    z.attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                    field_type
                );

                py::array_t<complex128> field_array = _vector_to_numpy(field);

                return field_array * ureg.attr("volt/meter");
            },
            py::arg("x"),
            py::arg("y"),
            py::arg("z"),
            py::arg("field_type"),
            R"pbdoc(
                Computes incident near-field electromagnetic fields.
                This method calculates the incident electromagnetic fields at specified points
                using the properties of the incident wave.
                Parameters
                ----------
                x : numpy.ndarray
                    Array of x coordinates of observation points.
                y : numpy.ndarray
                    Array of y coordinates of observation points.
                z : numpy.ndarray
                    Array of z coordinates of observation points.
                field_type : str
                    Field component type: "Ex", "Ey", "Ez", "\|E\|"
                source : BaseSource
                    The source of the electromagnetic field, used to determine incident field properties for near-field calculations.

                Returns
                -------
                numpy.ndarray
                    Complex array of field values at specified points.
                Notes
                -----
                This method does not depend on scatterer properties and computes fields
                based solely on the incident wave characteristics.
            )pbdoc"
        )

        .def(
            "get_total_nearfields",
            [ureg](
                Setup& self,
                const py::object& x,
                const py::object& y,
                const py::object& z,
                const std::string& field_type
            ) {
                std::vector<double> x_values = x.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();
                std::vector<double> y_values = y.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();
                std::vector<double> z_values = z.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                std::vector<complex128> field = self.get_total_nearfields(x_values, y_values, z_values, field_type);

                return vector_move_from_numpy(field, {field.size()}) * ureg.attr("volt/meter");
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
                    Field component type: "Ex", "Ey", "Ez", "Hx", "Hy", "Hz", "\|E\|", "\|H\|"

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
    .def(
        "get_representation",
        [](const Setup& self, const std::string& representation_type, py::kwargs kwargs) {
            std::string module_name;
            std::string class_name;

            if (representation_type == "farfields") {
                module_name = "PyMieSim.single.representations.farfields";
                class_name = "FarFields";
            }
            else if (representation_type == "stokes") {
                module_name = "PyMieSim.single.representations.stokes";
                class_name = "Stokes";
            }
            else if (representation_type == "spf") {
                module_name = "PyMieSim.single.representations.spf";
                class_name = "SPF";
            }
            else if (representation_type == "s1s2") {
                module_name = "PyMieSim.single.representations.s1s2";
                class_name = "S1S2";
            }
            else if (representation_type == "nearfields") {
                module_name = "PyMieSim.single.representations.nearfields";
                class_name = "NearFields";
            }
            else if (representation_type == "footprint") {
                module_name = "PyMieSim.single.representations.footprint";
                class_name = "Footprint";
            }
            else {
                throw std::runtime_error("Unknown representation type: " + representation_type);
            }

            py::module_ representation_module = py::module_::import(module_name.c_str());
            py::object representation_class = representation_module.attr(class_name.c_str());

            return representation_class(py::cast(self), **kwargs);
        },
        py::arg("representation_type")
    )
    .def(
        "plot_system",
        [](
            const Setup& self,
            const bool show_axes,
            const bool show_colorbar,
            const bool show_detector_cone,
            const bool show_unit_sphere,
            const double figure_size
        ) {
            py::module_ pyplot = py::module_::import("matplotlib.pyplot");

            py::object figure = pyplot.attr("figure")(
                py::arg("figsize") = py::make_tuple(figure_size, figure_size)
            );

            py::object ax = figure.attr("add_subplot")(
                py::int_(111),
                py::arg("projection") = py::str("3d")
            );

            py::cast(self.source).attr("_add_to_ax")(
                ax,
                py::arg("show_axes") = show_axes
            );

            self.scatterer->init(self.source);

            py::cast(self.scatterer).attr("_add_to_ax")(
                ax,
                py::arg("show_axes") = show_axes,
                py::arg("show_unit_sphere") = show_unit_sphere
            );

            if (self.detector) {
                self.detector->medium->initialize(self.source->wavelength);
                self.detector->initialize_mesh(self.scatterer);

                py::object detector_object = py::cast(self.detector);

                if (py::hasattr(detector_object, "_add_to_ax")) {
                    py::dict detector_kwargs;
                    detector_kwargs["show_axes"] = py::bool_(show_axes);

                    if (py::hasattr(detector_object, "mode_field")) {
                        detector_kwargs["show_colorbar"] = py::bool_(show_colorbar);
                        detector_kwargs["show_cone"] = py::bool_(show_detector_cone);
                    }
                    else {
                        detector_kwargs["show_cone"] = py::bool_(show_detector_cone);
                    }

                    detector_object.attr("_add_to_ax")(ax, **detector_kwargs);
                }
            }

            figure.attr("tight_layout")();
            pyplot.attr("show")();

            return figure;
        },
        py::arg("show_axes") = false,
        py::arg("show_colorbar") = true,
        py::arg("show_detector_cone") = false,
        py::arg("show_unit_sphere") = true,
        py::arg("figure_size") = 7.0,
        R"pbdoc(
            Plot the full single-scatterer optical system using Matplotlib.

            The plot includes the incident source, the scatterer, and the detector
            geometry when a detector is defined. The source is drawn using its
            propagation and polarization vectors. The scatterer is drawn using its
            normalized Matplotlib representation. The detector is initialized for
            the current scatterer and source before being drawn.

            Parameters
            ----------
            show_axes : bool, optional
                If ``True``, display Cartesian axis labels, ticks, and panes.
                If ``False``, hide the Matplotlib 3D axis frame after setting
                the plotting limits.
            show_colorbar : bool, optional
                If ``True``, show the detector colorbar when the detector supports
                field-colored plotting, such as ``CoherentMode``.
            show_detector_cone : bool, optional
                If ``True``, draw the detector collection cone when supported.
            show_unit_sphere : bool, optional
                If ``True``, draw the transparent unit sphere for scatterer
                angular reference when supported.
            figure_size : float, optional
                Width and height of the Matplotlib figure in inches.

            Returns
            -------
            matplotlib.figure.Figure
                Matplotlib figure containing the system visualization.

            Notes
            -----
            This method replaces the previous PyVista-based visualization path.
            It uses the private ``_add_to_ax`` helpers exposed by source,
            scatterer, and detector bindings.
        )pbdoc"
    )
;


}
