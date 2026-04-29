#include <algorithm>
#include <cmath>
#include <complex>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <pint/pint.h>

#include <utils/casting.h>
#include <utils/numpy_interface.h>

#include <single/scatterer/utils.h>
#include <single/detector/photodiode.h>
#include <single/detector/coherent_mode.h>
#include <single/detector/integrating_sphere.h>


namespace py = pybind11;


PYBIND11_MODULE(detector, module) {
    py::object ureg = get_shared_ureg();

    module.doc() = R"pbdoc(
        Detector bindings for PyMieSim.

        This module exposes detector classes used to compute optical coupling
        between Lorenz Mie scatterers and detector geometries. It includes
        photodiode, coherent mode, and integrating sphere detectors.

        The plotting helpers use Matplotlib axes directly through private
        ``_add_to_ax`` methods. Figure creation, styling, and public plotting
        APIs should remain controlled from the Python layer.
    )pbdoc";

    auto ensure_matplotlib_3d_axis = [](
        const py::object& ax
    ) -> void {
        if (!py::hasattr(ax, "get_zlim")) {
            throw std::runtime_error(
                "Detector plotting requires a Matplotlib 3D axis. "
                "Create it with: figure.add_subplot(111, projection=\"3d\")."
            );
        }
    };

    auto numpy_array_from_vector = [](
        const std::vector<double>& values
    ) -> py::array_t<double> {
        py::array_t<double> array(
            static_cast<py::ssize_t>(values.size())
        );

        auto mutable_array = array.mutable_unchecked<1>();

        for (py::ssize_t index = 0; index < static_cast<py::ssize_t>(values.size()); ++index) {
            mutable_array(index) = values[static_cast<std::size_t>(index)];
        }

        return array;
    };

    auto numpy_surface_from_vector = [](
        const std::vector<double>& values,
        const std::size_t radial_sampling,
        const std::size_t angular_sampling
    ) -> py::array_t<double> {
        py::array_t<double> array(
            py::array::ShapeContainer{
                static_cast<py::ssize_t>(radial_sampling),
                static_cast<py::ssize_t>(angular_sampling)
            }
        );

        auto mutable_array = array.mutable_unchecked<2>();

        for (std::size_t radial_index = 0; radial_index < radial_sampling; ++radial_index) {
            for (std::size_t angular_index = 0; angular_index < angular_sampling; ++angular_index) {
                const std::size_t flat_index = radial_index * angular_sampling + angular_index;

                mutable_array(
                    static_cast<py::ssize_t>(radial_index),
                    static_cast<py::ssize_t>(angular_index)
                ) = values[flat_index];
            }
        }

        return array;
    };

    auto set_detector_axis_limits = [](
        const py::object& ax,
        const bool show_axes
    ) -> void {
        ax.attr("set_xlim")(-1.0, 1.0);
        ax.attr("set_ylim")(-1.0, 1.0);
        ax.attr("set_zlim")(-1.0, 1.0);

        if (py::hasattr(ax, "set_box_aspect")) {
            ax.attr("set_box_aspect")(py::make_tuple(1.0, 1.0, 1.0));
        }

        if (show_axes) {
            ax.attr("set_xlabel")("x");
            ax.attr("set_ylabel")("y");
            ax.attr("set_zlabel")("z");
        } else {
            ax.attr("set_axis_off")();
        }
    };

    auto add_collection_cone_to_ax = [
        numpy_surface_from_vector
    ](
        const BaseDetector& detector,
        const py::object& ax,
        const py::object& cone_color,
        const double cone_alpha,
        const std::size_t cone_angular_sampling,
        const std::size_t cone_radial_sampling
    ) -> void {
        const BaseDetector::CollectionConeSurface cone_surface =
            detector.get_collection_cone_surface(
                cone_angular_sampling,
                cone_radial_sampling
            );

        const std::size_t angular_sampling = cone_surface.angular_sampling;
        const std::size_t rim_offset =
            (cone_surface.radial_sampling - 1) * angular_sampling;

        std::vector<double> lateral_x(2 * angular_sampling, 0.0);
        std::vector<double> lateral_y(2 * angular_sampling, 0.0);
        std::vector<double> lateral_z(2 * angular_sampling, 0.0);

        for (std::size_t angular_index = 0; angular_index < angular_sampling; ++angular_index) {
            const std::size_t rim_index = rim_offset + angular_index;
            const std::size_t lateral_index = angular_sampling + angular_index;

            lateral_x[lateral_index] = cone_surface.x[rim_index];
            lateral_y[lateral_index] = cone_surface.y[rim_index];
            lateral_z[lateral_index] = cone_surface.z[rim_index];
        }

        ax.attr("plot_surface")(
            numpy_surface_from_vector(lateral_x, 2, angular_sampling),
            numpy_surface_from_vector(lateral_y, 2, angular_sampling),
            numpy_surface_from_vector(lateral_z, 2, angular_sampling),
            py::arg("color") = cone_color,
            py::arg("alpha") = cone_alpha,
            py::arg("linewidth") = 0.0,
            py::arg("edgecolor") = py::str("none"),
            py::arg("shade") = false,
            py::arg("antialiased") = false,
            py::arg("zorder") = 0
        );
    };

    auto add_detector_points_to_ax = [
        numpy_array_from_vector
    ](
        const BaseDetector& detector,
        const py::object& ax,
        const py::object& color,
        const double field_point_size,
        const double point_radial_offset
    ) -> py::object {
        const BaseDetector::CartesianCoordinateVectors coordinates =
            detector.get_cartesian_coordinate_vectors();

        std::vector<double> x_values = coordinates.x;
        std::vector<double> y_values = coordinates.y;
        std::vector<double> z_values = coordinates.z;

        for (std::size_t index = 0; index < x_values.size(); ++index) {
            x_values[index] *= point_radial_offset;
            y_values[index] *= point_radial_offset;
            z_values[index] *= point_radial_offset;
        }

        py::array_t<double> x_array = numpy_array_from_vector(x_values);
        py::array_t<double> y_array = numpy_array_from_vector(y_values);
        py::array_t<double> z_array = numpy_array_from_vector(z_values);

        return ax.attr("scatter")(
            x_array,
            y_array,
            py::arg("zs") = z_array,
            py::arg("c") = color,
            py::arg("s") = field_point_size,
            py::arg("depthshade") = false,
            py::arg("zorder") = 10
        );
    };

    auto add_colored_detector_points_to_ax = [
        numpy_array_from_vector
    ](
        const BaseDetector& detector,
        const py::object& ax,
        const std::vector<double>& scalar_values,
        const double field_point_size,
        const double point_radial_offset,
        const py::object& colormap,
        const double color_minimum,
        const double color_maximum
    ) -> py::object {
        const BaseDetector::CartesianCoordinateVectors coordinates =
            detector.get_cartesian_coordinate_vectors();

        std::vector<double> x_values = coordinates.x;
        std::vector<double> y_values = coordinates.y;
        std::vector<double> z_values = coordinates.z;

        for (std::size_t index = 0; index < x_values.size(); ++index) {
            x_values[index] *= point_radial_offset;
            y_values[index] *= point_radial_offset;
            z_values[index] *= point_radial_offset;
        }

        py::array_t<double> x_array = numpy_array_from_vector(x_values);
        py::array_t<double> y_array = numpy_array_from_vector(y_values);
        py::array_t<double> z_array = numpy_array_from_vector(z_values);
        py::array_t<double> scalar_array = numpy_array_from_vector(scalar_values);

        return ax.attr("scatter")(
            x_array,
            y_array,
            py::arg("zs") = z_array,
            py::arg("c") = scalar_array,
            py::arg("s") = field_point_size,
            py::arg("cmap") = colormap,
            py::arg("vmin") = color_minimum,
            py::arg("vmax") = color_maximum,
            py::arg("depthshade") = false,
            py::arg("zorder") = 10
        );
    };

    py::class_<BaseDetector, std::shared_ptr<BaseDetector>>(module, "BaseDetector")
        .def_readwrite(
            "medium",
            &BaseDetector::medium,
            R"pbdoc(
                Medium in which the detector is immersed.

                The medium determines the refractive index surrounding the
                detector and is used when computing angular collection and
                interface transmission effects.
            )pbdoc"
        )
        .def_property_readonly(
            "numerical_aperture",
            [](
                BaseDetector& self
            ) {
                return py::float_(self.numerical_aperture);
            },
            R"pbdoc(
                Numerical aperture of the detector.

                Controls the angular extent of the collected scattered field.
            )pbdoc"
        )
        .def_readonly(
            "mesh",
            &BaseDetector::fibonacci_mesh,
            R"pbdoc(
                Fibonacci angular mesh used by the detector.

                The mesh defines the angular directions on which detector fields
                and coupling integrals are sampled.
            )pbdoc"
        )
        .def_property_readonly(
            "max_angle",
            [ureg](
                BaseDetector& self
            ) {
                return py::float_(self.max_angle) * ureg.attr("radian");
            },
            R"pbdoc(
                Maximum detector collection angle.

                Returns
                -------
                pint.Quantity
                    Maximum collection angle in radians.
            )pbdoc"
        )
        .def_property_readonly(
            "min_angle",
            [ureg](
                BaseDetector& self
            ) {
                return py::float_(self.min_angle) * ureg.attr("radian");
            },
            R"pbdoc(
                Minimum detector collection angle.

                Returns
                -------
                pint.Quantity
                    Minimum collection angle in radians.
            )pbdoc"
        )
        .def(
            "initialize_mesh",
            &BaseDetector::initialize_mesh,
            py::arg("scatterer"),
            R"pbdoc(
                Initialize the detector mesh for a scatterer.

                Parameters
                ----------
                scatterer : BaseScatterer
                    Scatterer used to configure detector dependent angular
                    quantities.
            )pbdoc"
        )
        .def_property(
            "scalar_field",
            [](
                BaseDetector& self
            ) {
                return vector_as_numpy_view(self, self.scalar_field);
            },
            [](
                BaseDetector& self,
                py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> array
            ) {
                vector_assign_from_numpy(self.scalar_field, array);
            },
            R"pbdoc(
                Complex detector scalar field samples.

                The getter returns a zero copy NumPy view tied to the detector
                lifetime. The setter accepts any one dimensional array castable to
                ``numpy.complex128``.

                Returns
                -------
                numpy.ndarray
                    Complex array with shape ``(sampling,)``.
            )pbdoc"
        )
        .def(
            "get_structured_scalarfield",
            [](
                BaseDetector& self,
                const size_t sampling
            ) {
                std::vector<complex128> vector = self.get_structured_scalarfield(sampling);
                return vector_move_from_numpy(std::move(vector), {sampling, sampling});
            },
            py::arg("sampling"),
            R"pbdoc(
                Compute a structured representation of the scalar field.

                Parameters
                ----------
                sampling : int
                    Number of samples per structured angular dimension.

                Returns
                -------
                numpy.ndarray
                    Complex scalar field array with shape ``(sampling, sampling)``.
            )pbdoc"
        )
        .def(
            "get_poynting_field",
            [ureg](
                BaseDetector& self,
                std::shared_ptr<BaseScatterer> scatterer,
                const py::object& distance,
                std::shared_ptr<BaseSource> source
            ) {
                const double distance_value = distance.attr("to")("meter").attr("magnitude").cast<double>();

                std::vector<double> vector = self.get_poynting_field(
                    scatterer,
                    source,
                    distance_value
                );

                py::array_t<double> array = vector_move_from_numpy(
                    std::move(vector),
                    {vector.size()}
                );

                return (array * ureg.attr("watt/meter**2")).attr("to_compact")();
            },
            py::arg("scatterer"),
            py::arg("distance") = 1.0,
            py::arg("source"),
            R"pbdoc(
                Compute the sampled Poynting field magnitude.

                Parameters
                ----------
                scatterer : BaseScatterer
                    Scatterer that produces the scattered field.
                distance : pint.Quantity
                    Observation distance from the scatterer.
                source : BaseSource
                    Incident source illuminating the scatterer.

                Returns
                -------
                pint.Quantity
                    Poynting field magnitude sampled on the detector mesh, with
                    units of watts per square meter.
            )pbdoc"
        )
        .def(
            "get_energy_flow",
            [ureg](
                BaseDetector& self,
                std::shared_ptr<BaseScatterer> scatterer,
                std::shared_ptr<BaseSource> source
            ) {
                const double energy_flow = self.get_energy_flow(scatterer, source);
                return (py::float_(energy_flow) * ureg.attr("watt")).attr("to_compact")();
            },
            py::arg("scatterer"),
            py::arg("source"),
            R"pbdoc(
                Compute the total radiated energy flow sampled by the detector.

                Parameters
                ----------
                scatterer : BaseScatterer
                    Scatterer that produces the scattered field.
                source : BaseSource
                    Incident source illuminating the scatterer.

                Returns
                -------
                pint.Quantity
                    Total energy flow in watts.
            )pbdoc"
        )
        .def_readonly(
            "sampling",
            &BaseDetector::sampling,
            R"pbdoc(
                Number of angular samples used by the detector mesh.
            )pbdoc"
        )
        .def_property_readonly(
            "phi_offset",
            [ureg](
                BaseDetector& self
            ) {
                return py::float_(self.phi_offset) * ureg.attr("radian");
            },
            R"pbdoc(
                Azimuthal orientation offset of the detector.

                Returns
                -------
                pint.Quantity
                    Azimuthal offset in radians.
            )pbdoc"
        )
        .def_property_readonly(
            "gamma_offset",
            [ureg](
                BaseDetector& self
            ) {
                return py::float_(self.gamma_offset) * ureg.attr("radian");
            },
            R"pbdoc(
                Polar orientation offset of the detector.

                Returns
                -------
                pint.Quantity
                    Polar offset in radians.
            )pbdoc"
        )
        .def_property_readonly(
            "polarization_filter",
            [](
                BaseDetector& self
            ) {
                return py::cast(self.polarization_filter);
            },
            R"pbdoc(
                Polarization filter applied before coupling evaluation.
            )pbdoc"
        )
        .def_readonly(
            "mode_field",
            &BaseDetector::mode_field,
            R"pbdoc(
                Mode field associated with coherent detector calculations.
            )pbdoc"
        )
        .def_readonly(
            "mode_id",
            &BaseDetector::mode_id,
            R"pbdoc(
                Identifier of the detector mode family and order.
            )pbdoc"
        )
        .def(
            "get_coupling",
            [ureg](
                BaseDetector& self,
                std::shared_ptr<BaseScatterer> scatterer,
                std::shared_ptr<BaseSource> source
            ) {
                const double coupling = self.get_coupling(scatterer, source);
                return (py::float_(coupling) * ureg.attr("watt")).attr("to_compact")();
            },
            py::arg("scatterer"),
            py::arg("source"),
            R"pbdoc(
                Compute optical power coupled into the detector.

                Parameters
                ----------
                scatterer : BaseScatterer
                    Scatterer producing the scattered field.
                source : BaseSource
                    Incident source illuminating the scatterer.

                Returns
                -------
                pint.Quantity
                    Coupled optical power in watts.
            )pbdoc"
        );

    py::class_<Photodiode, BaseDetector, std::shared_ptr<Photodiode>>(
        module,
        "Photodiode",
        R"pbdoc(
            Photodiode detector for Lorenz Mie scattering simulations.

            A photodiode integrates the scattered power over an angular
            collection aperture defined by its numerical aperture, orientation,
            optional obscuration cache, and polarization filter.
        )pbdoc"
    )
        .def(
            py::init(
                [ureg](
                    const py::object& numerical_aperture,
                    const py::object& phi_offset,
                    const py::object& gamma_offset,
                    const py::object& medium,
                    const py::object& polarization_filter,
                    const py::object& cache_numerical_aperture,
                    const std::size_t sampling
                ) {
                    const double numerical_aperture_value = numerical_aperture.cast<double>();

                    std::shared_ptr<BaseMedium> parsed_medium;

                    if (medium.is(py::none())) {
                        parsed_medium = std::make_shared<ConstantMedium>(1.0);
                    } else {
                        parsed_medium = parse_medium_object(medium, ureg);
                    }

                    return std::make_shared<Photodiode>(
                        sampling,
                        numerical_aperture_value,
                        cache_numerical_aperture.cast<double>(),
                        phi_offset.attr("to")("radian").attr("magnitude").cast<double>(),
                        gamma_offset.attr("to")("radian").attr("magnitude").cast<double>(),
                        Casting::Polarization::cast_py_to_polarization_state(polarization_filter),
                        std::move(parsed_medium)
                    );
                }
            ),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("medium") = py::none(),
            py::arg("polarization_filter") = PolarizationState(),
            py::arg("cache_numerical_aperture") = py::float_(0.0),
            py::arg("sampling") = 200,
            R"pbdoc(
                Construct a photodiode detector.

                Parameters
                ----------
                numerical_aperture : float
                    Numerical aperture defining the detector collection angle.
                phi_offset : pint.Quantity
                    Azimuthal orientation offset.
                gamma_offset : pint.Quantity
                    Polar orientation offset.
                medium : BaseMedium or float, optional
                    Surrounding medium. If ``None``, a constant refractive index
                    of 1.0 is used.
                polarization_filter : PolarizationState, optional
                    Polarization filter applied before detection.
                cache_numerical_aperture : float, optional
                    Numerical aperture of the central obscuration cache.
                sampling : int, optional
                    Number of angular samples on the detector mesh.
            )pbdoc"
        )
        .def(
            "print_properties",
            &Photodiode::print_properties,
            py::arg("precision") = 4,
            R"pbdoc(
                Print photodiode properties.

                Parameters
                ----------
                precision : int, optional
                    Number of decimal places used for floating point values.
            )pbdoc"
        )
        .def(
            "_add_to_ax",
            [
                ensure_matplotlib_3d_axis,
                add_detector_points_to_ax,
                add_collection_cone_to_ax,
                set_detector_axis_limits
            ](
                const Photodiode& self,
                const py::object& ax,
                const py::object& cone_color,
                const double field_point_size,
                const double point_radial_offset,
                const double cone_alpha,
                const bool show_cone,
                const bool show_axes,
                const std::size_t cone_angular_sampling,
                const std::size_t cone_radial_sampling
            ) -> void {
                ensure_matplotlib_3d_axis(ax);

                if (show_cone) {
                    add_collection_cone_to_ax(
                        self,
                        ax,
                        cone_color,
                        cone_alpha,
                        cone_angular_sampling,
                        cone_radial_sampling
                    );
                }

                add_detector_points_to_ax(
                    self,
                    ax,
                    py::str("black"),
                    field_point_size,
                    point_radial_offset
                );

                set_detector_axis_limits(ax, show_axes);
            },
            py::arg("ax"),
            py::arg("cone_color") = py::str("red"),
            py::arg("field_point_size") = 20.0,
            py::arg("point_radial_offset") = 1.025,
            py::arg("cone_alpha") = 0.20,
            py::arg("show_cone") = true,
            py::arg("show_axes") = false,
            py::arg("cone_angular_sampling") = 24,
            py::arg("cone_radial_sampling") = 12,
            R"pbdoc(
                Add the photodiode detector geometry to a Matplotlib 3D axis.

                This private helper draws the detector angular sampling directions
                as points on the unit sphere and, optionally, draws a translucent
                lateral collection cone representing the angular aperture. The
                circular base of the cone is not drawn.

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Matplotlib axis created with ``projection="3d"``.
                cone_color : str or tuple, optional
                    Matplotlib color used for the collection cone.
                field_point_size : float, optional
                    Marker size used for detector sampling points.
                point_radial_offset : float, optional
                    Multiplicative radial offset applied to detector sampling
                    points for visualization. Values slightly above 1 place the
                    points just outside the collection cone surface.
                cone_alpha : float, optional
                    Transparency of the collection cone.
                show_cone : bool, optional
                    If ``True``, draw the collection cone.
                show_axes : bool, optional
                    If ``True``, display Cartesian axis labels, ticks, and panes.
                    If ``False``, hide the Matplotlib 3D axis frame after setting
                    the plotting limits.
                cone_angular_sampling : int, optional
                    Number of angular samples used to draw the cone surface.
                cone_radial_sampling : int, optional
                    Number of radial samples used to generate the cone geometry.

                Notes
                -----
                This method replaces the previous PyVista based ``add_to_scene``
                method. It operates on an existing Matplotlib axis so the Python
                plotting layer remains responsible for figure creation and style.
            )pbdoc"
        );

    py::class_<CoherentMode, BaseDetector, std::shared_ptr<CoherentMode>>(
        module,
        "CoherentMode",
        R"pbdoc(
            Coherent mode detector for Lorenz Mie scattering simulations.

            A coherent mode detector projects the scattered field onto a specified
            collecting mode and can account for orientation, rotation, numerical
            aperture, obscuration cache, and polarization filtering.
        )pbdoc"
    )
        .def(
            py::init(
                [ureg](
                    const std::string& mode_number,
                    const py::object& numerical_aperture,
                    const py::object& phi_offset,
                    const py::object& gamma_offset,
                    const py::object& rotation,
                    const py::object& medium,
                    const py::object& cache_numerical_aperture,
                    const py::object& polarization_filter,
                    const py::object& mean_coupling,
                    const std::size_t& sampling
                ) {
                    const double numerical_aperture_value = numerical_aperture.cast<double>();

                    std::shared_ptr<BaseMedium> parsed_medium;

                    if (medium.is(py::none())) {
                        parsed_medium = std::make_shared<ConstantMedium>(1.0);
                    } else {
                        parsed_medium = parse_medium_object(medium, ureg);
                    }

                    if (numerical_aperture_value < 0.0) {
                        throw std::runtime_error("Numerical aperture must be non negative.");
                    }

                    if (numerical_aperture_value > 0.342) {
                        std::printf(
                            "Warning: numerical aperture exceeds 0.342 for coherent detectors, which may violate the paraxial approximation.\n"
                        );
                    }

                    return std::make_shared<CoherentMode>(
                        mode_number,
                        sampling,
                        numerical_aperture_value,
                        cache_numerical_aperture.cast<double>(),
                        phi_offset.attr("to")("radian").attr("magnitude").cast<double>(),
                        gamma_offset.attr("to")("radian").attr("magnitude").cast<double>(),
                        Casting::Polarization::cast_py_to_polarization_state(polarization_filter),
                        rotation.attr("to")("radian").attr("magnitude").cast<double>(),
                        mean_coupling.cast<bool>(),
                        std::move(parsed_medium)
                    );
                }
            ),
            py::arg("mode_number"),
            py::arg("numerical_aperture"),
            py::arg("phi_offset"),
            py::arg("gamma_offset"),
            py::arg("rotation"),
            py::arg("medium") = py::none(),
            py::arg("cache_numerical_aperture") = py::float_(0.0),
            py::arg("polarization_filter") = PolarizationState(),
            py::arg("mean_coupling") = false,
            py::arg("sampling") = 200,
            R"pbdoc(
                Construct a coherent mode detector.

                Parameters
                ----------
                mode_number : str
                    Mode identifier, for example ``"LP01"`` or ``"HG12"``.
                numerical_aperture : float
                    Numerical aperture defining the angular collection region.
                phi_offset : pint.Quantity
                    Azimuthal orientation offset.
                gamma_offset : pint.Quantity
                    Polar orientation offset.
                rotation : pint.Quantity
                    Rotation of the detector mode field.
                medium : BaseMedium or float, optional
                    Surrounding medium. If ``None``, a constant refractive index
                    of 1.0 is used.
                cache_numerical_aperture : float, optional
                    Numerical aperture of the central obscuration cache.
                polarization_filter : PolarizationState, optional
                    Polarization filter applied before coherent projection.
                mean_coupling : bool, optional
                    If ``True``, use mean coupling behavior where supported by
                    the detector implementation.
                sampling : int, optional
                    Number of angular samples on the detector mesh.
            )pbdoc"
        )
        .def_property_readonly(
            "rotation",
            [ureg](
                CoherentMode& self
            ) {
                return py::float_(self.rotation) * ureg.attr("radian");
            },
            R"pbdoc(
                Rotation angle of the coherent mode field.

                Returns
                -------
                pint.Quantity
                    Rotation angle in radians.
            )pbdoc"
        )
        .def(
            "print_properties",
            &CoherentMode::print_properties,
            py::arg("precision") = 4,
            R"pbdoc(
                Print coherent mode detector properties.

                Parameters
                ----------
                precision : int, optional
                    Number of decimal places used for floating point values.
            )pbdoc"
        )
        .def(
            "_add_to_ax",
            [
                ensure_matplotlib_3d_axis,
                add_colored_detector_points_to_ax,
                add_collection_cone_to_ax,
                set_detector_axis_limits
            ](
                const CoherentMode& self,
                const py::object& ax,
                const py::object& cone_color,
                const double field_point_size,
                const double point_radial_offset,
                const double cone_alpha,
                const bool show_cone,
                const bool show_colorbar,
                const bool show_axes,
                const std::size_t cone_angular_sampling,
                const std::size_t cone_radial_sampling
            ) -> void {
                ensure_matplotlib_3d_axis(ax);

                py::object blue_black_red =
                    py::module_::import("MPSPlots.colormaps").attr("blue_black_red");

                std::vector<double> scalar_field_real_values;
                scalar_field_real_values.resize(self.sampling, 0.0);

                double absolute_maximum = 0.0;

                const std::size_t common_size = std::min(
                    self.scalar_field.size(),
                    static_cast<std::size_t>(self.sampling)
                );

                for (std::size_t index = 0; index < common_size; ++index) {
                    const double value = self.scalar_field[index].real();

                    scalar_field_real_values[index] = value;
                    absolute_maximum = std::max(absolute_maximum, std::abs(value));
                }

                if (absolute_maximum <= 0.0) {
                    absolute_maximum = 1.0;
                }

                if (show_cone) {
                    add_collection_cone_to_ax(
                        self,
                        ax,
                        cone_color,
                        cone_alpha,
                        cone_angular_sampling,
                        cone_radial_sampling
                    );
                }

                py::object scatter = add_colored_detector_points_to_ax(
                    self,
                    ax,
                    scalar_field_real_values,
                    field_point_size,
                    point_radial_offset,
                    blue_black_red,
                    -absolute_maximum,
                    absolute_maximum
                );

                if (show_colorbar) {
                    ax.attr("figure").attr("colorbar")(
                        scatter,
                        py::arg("ax") = ax,
                        py::arg("label") = "Collecting field real part"
                    );
                }

                set_detector_axis_limits(ax, show_axes);
            },
            py::arg("ax"),
            py::arg("cone_color") = py::str("red"),
            py::arg("field_point_size") = 20.0,
            py::arg("point_radial_offset") = 1.025,
            py::arg("cone_alpha") = 0.20,
            py::arg("show_cone") = true,
            py::arg("show_colorbar") = true,
            py::arg("show_axes") = false,
            py::arg("cone_angular_sampling") = 96,
            py::arg("cone_radial_sampling") = 16,
            R"pbdoc(
                Add the coherent detector mode to a Matplotlib 3D axis.

                This private helper draws the detector angular sampling directions
                on the unit sphere and colors each point by the real part of the
                coherent collecting field. The color scale is centered on zero so
                positive and negative field regions are displayed symmetrically.
                The collection cone is drawn before the points, and the points are
                radially offset for visualization so they remain visible above the
                cone surface.

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Matplotlib axis created with ``projection="3d"``.
                cone_color : str or tuple, optional
                    Matplotlib color used for the collection cone.
                field_point_size : float, optional
                    Marker size used for detector sampling points.
                point_radial_offset : float, optional
                    Multiplicative radial offset applied to detector sampling
                    points for visualization. Values slightly above 1 place the
                    points just outside the collection cone surface.
                cone_alpha : float, optional
                    Transparency of the collection cone.
                show_cone : bool, optional
                    If ``True``, draw the lateral collection cone.
                show_colorbar : bool, optional
                    If ``True``, add a colorbar for the real part of the
                    collecting field.
                show_axes : bool, optional
                    If ``True``, display Cartesian axis labels, ticks, and panes.
                    If ``False``, hide the Matplotlib 3D axis frame after setting
                    the plotting limits.
                cone_angular_sampling : int, optional
                    Number of angular samples used to draw the cone surface.
                cone_radial_sampling : int, optional
                    Number of radial samples used to generate the cone geometry.

                Notes
                -----
                This method replaces the previous PyVista based ``add_to_scene``
                method. It uses ``MPSPlots.colormaps.blue_black_red`` as the
                diverging colormap for the real part of the coherent collecting
                field.
            )pbdoc"
        );

    py::class_<IntegratingSphere, BaseDetector, std::shared_ptr<IntegratingSphere>>(
        module,
        "IntegratingSphere",
        R"pbdoc(
            Integrating sphere detector for Lorenz Mie scattering simulations.

            The integrating sphere represents full angular collection over the
            detector mesh, optionally with polarization filtering.
        )pbdoc"
    )
        .def(
            py::init(
                [](
                    const py::object& polarization_filter,
                    const std::size_t& sampling
                ) {
                    return std::make_shared<IntegratingSphere>(
                        sampling,
                        Casting::Polarization::cast_py_to_polarization_state(polarization_filter)
                    );
                }
            ),
            py::arg("polarization_filter") = PolarizationState(),
            py::arg("sampling") = 400,
            R"pbdoc(
                Construct an integrating sphere detector.

                Parameters
                ----------
                polarization_filter : PolarizationState, optional
                    Polarization filter applied before detection.
                sampling : int, optional
                    Number of angular samples on the detector mesh.
            )pbdoc"
        )
        .def(
            "print_properties",
            &IntegratingSphere::print_properties,
            py::arg("precision") = 4,
            R"pbdoc(
                Print integrating sphere detector properties.

                Parameters
                ----------
                precision : int, optional
                    Number of decimal places used for floating point values.
            )pbdoc"
        )
        .def(
            "_add_to_ax",
            [
                ensure_matplotlib_3d_axis,
                add_detector_points_to_ax,
                set_detector_axis_limits
            ](
                const IntegratingSphere& self,
                const py::object& ax,
                const double field_point_size,
                const double point_radial_offset,
                const bool show_axes
            ) -> void {
                ensure_matplotlib_3d_axis(ax);

                add_detector_points_to_ax(
                    self,
                    ax,
                    py::str("black"),
                    field_point_size,
                    point_radial_offset
                );

                set_detector_axis_limits(ax, show_axes);
            },
            py::arg("ax"),
            py::arg("field_point_size") = 20.0,
            py::arg("point_radial_offset") = 1.025,
            py::arg("show_axes") = false,
            R"pbdoc(
                Add the integrating sphere angular mesh to a Matplotlib 3D axis.

                This private helper draws the angular sampling directions of the
                integrating sphere as points on the unit sphere. No collection
                cone is drawn because the integrating sphere represents full
                angular collection.

                Parameters
                ----------
                ax : matplotlib.axes.Axes
                    Matplotlib axis created with ``projection="3d"``.
                field_point_size : float, optional
                    Marker size used for detector sampling points.
                point_radial_offset : float, optional
                    Multiplicative radial offset applied to detector sampling
                    points for visualization.
                show_axes : bool, optional
                    If ``True``, display Cartesian axis labels, ticks, and panes.
                    If ``False``, hide the Matplotlib 3D axis frame after setting
                    the plotting limits.

                Notes
                -----
                This method replaces the previous PyVista based ``add_to_scene``
                method. It operates on an existing Matplotlib axis so the Python
                plotting layer remains responsible for figure creation and style.
            )pbdoc"
        );
}