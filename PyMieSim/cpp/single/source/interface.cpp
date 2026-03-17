#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include <array>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>


#include <pint/pint.h>
#include <utils/numpy_interface.h>
#include "./source.h"
#include <single/utils.h>
#include <utils/casting.h>

namespace py = pybind11;

PYBIND11_MODULE(source, module)
{
    py::object ureg = get_shared_ureg();
    module.doc() = "Source bindings for PyMieSim.";

    py::class_<BaseSource, std::shared_ptr<BaseSource>>(module, "BaseSource",
        R"pbdoc(
            Base class for optical excitation sources.

            A ``BaseSource`` defines the incident field that illuminates a scatterer.
            Concrete sources provide the wavelength, electric field amplitude, and
            polarization state.

            Notes
            -----
            Units are returned as Pint quantities using the shared unit registry.
            Polarization is represented internally as a Jones vector with shape ``(N, 2)``,
            stored as complex values.
        )pbdoc"
        )
        .def(py::init<>(),
            R"pbdoc(
                Construct a source base instance.

                Notes
                -----
                This constructor is primarily used for binding completeness.
                In user code, you typically instantiate a concrete source such as
                :class:`~PyMieSim.single.source.PlaneWave` or
                :class:`~PyMieSim.single.source.Gaussian`.
            )pbdoc"
        )
        .def_property_readonly(
            "wavelength",
            [ureg](const BaseSource& self) {
                return (py::float_(self.wavelength) * ureg.attr("meter"));
            },
            R"pbdoc(
                Wavelength in vacuum.

                Returns
                -------
                pint.Quantity
                    Wavelength with units of length.
            )pbdoc"
        )
        .def_property_readonly(
            "wavenumber_vacuum",
            [ureg](const BaseSource& self) {
                return (py::float_(self.wavenumber_vacuum) * ureg.attr("1/meter")).attr("to_compact")();
            },
            R"pbdoc(
                Vacuum wavenumber.

                This is defined as ``k0 = 2*pi / wavelength``.

                Returns
                -------
                pint.Quantity
                    Wavenumber with units of inverse length.
            )pbdoc"
        )
        .def_property_readonly(
            "amplitude",
            [ureg](const BaseSource& self) {
                return (py::float_(self.amplitude) * ureg.attr("volt/meter")).attr("to_compact")();
            },
            R"pbdoc(
                Electric field amplitude.

                Returns
                -------
                pint.Quantity
                    Field amplitude with units of electric field (V/m).
            )pbdoc"
        )
        .def_readonly(
            "polarization",
            &BaseSource::polarization,
            R"pbdoc(
                Convenience Jones vector.

                Returns the first row of :attr:`polarization_element` as a vector
                ``[Ex, Ey]``.

                Returns
                -------
                numpy.ndarray
                    Complex array of shape ``(2,)``.
            )pbdoc"
        )
        .def(
            "add_to_scene",
            [](const BaseSource& self, const py::object& scene, py::kwargs kwargs) {
                py::module_ pyvista = py::module_::import("pyvista");
                py::module_ numpy = py::module_::import("numpy");

                py::object Arrow = pyvista.attr("Arrow");

                const std::complex<double> Ex = self.polarization.jones_vector[0];
                const std::complex<double> Ey = self.polarization.jones_vector[1];

                const double electric_x = std::real(Ex);
                const double electric_y = std::real(Ey);

                const double electric_norm = std::sqrt(
                    electric_x * electric_x +
                    electric_y * electric_y
                );

                const double normalized_electric_x = (electric_norm > 0.0) ? electric_x / electric_norm : 1.0;
                const double normalized_electric_y = (electric_norm > 0.0) ? electric_y / electric_norm : 0.0;
                const double normalized_electric_z = 0.0;

                py::object k_arrow = Arrow(
                    py::arg("start") = py::make_tuple(0.0, 0.0, -2.1),
                    py::arg("direction") = py::make_tuple(0.0, 0.0, 1.0),
                    py::arg("tip_length") = py::float_(0.2),
                    py::arg("tip_radius") = py::float_(0.05),
                    py::arg("shaft_radius") = py::float_(0.02)
                );

                py::object electric_arrow = Arrow(
                    py::arg("start") = py::make_tuple(0.0, 0.0, -2.1),
                    py::arg("direction") = py::make_tuple(
                        normalized_electric_x,
                        normalized_electric_y,
                        normalized_electric_z
                    ),
                    py::arg("scale") = py::float_(0.8),
                    py::arg("tip_length") = py::float_(0.25),
                    py::arg("tip_radius") = py::float_(0.06),
                    py::arg("shaft_radius") = py::float_(0.025)
                );

                py::dict k_kwargs;
                k_kwargs["color"] = py::str("orange");
                if (kwargs.contains("opacity")) {
                    k_kwargs["opacity"] = kwargs["opacity"];
                }

                py::dict e_kwargs;
                e_kwargs["color"] = py::str("deepskyblue");
                if (kwargs.contains("opacity")) {
                    e_kwargs["opacity"] = kwargs["opacity"];
                }

                scene.attr("add_mesh")(k_arrow, **k_kwargs);
                scene.attr("add_mesh")(electric_arrow, **e_kwargs);
            },
            py::arg("scene")
        )
        ;


    py::class_<Planewave, BaseSource, std::shared_ptr<Planewave>>(module, "PlaneWave",
        R"pbdoc(
            Monochromatic plane wave source.

            This source models a spatially uniform incident field with a specified
            wavelength, polarization, and electric field amplitude.
        )pbdoc")
        .def(
            py::init([](
                const py::object& wavelength,
                const py::object& polarization,
                const py::object& amplitude
            ) {
                return std::make_shared<Planewave>(
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>(),
                    Casting::cast_py_to_polarization_state(polarization),
                    amplitude.attr("to")("volt/meter").attr("magnitude").cast<double>()
                );
            }),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("amplitude"),
            R"pbdoc(
                Construct a plane wave source.

                Parameters
                ----------
                wavelength : pint.Quantity
                    Wavelength in vacuum. Must be convertible to meters.
                polarization : object
                    Polarization specification. Accepted inputs are:

                    * A bound C++ ``JonesVector`` instance.
                    * A ``PyMieSim.single.polarization.BasePolarization`` instance.
                    * Any object accepted by ``PyMieSim.single.polarization.Linear(...)``
                    (it will be converted to a linear polarization).

                    The polarization is internally stored as a Jones matrix with shape
                    ``(N, 2)`` of complex values.
                amplitude : pint.Quantity
                    Electric field amplitude. Must be convertible to V/m.

                Returns
                -------
                PlaneWave
                    Configured plane wave source.

                Notes
                -----
                The incident field phase convention follows the C++ implementation.
            )pbdoc"
        );


    py::class_<Gaussian, BaseSource, std::shared_ptr<Gaussian>>(module, "Gaussian",
        R"pbdoc(
            Monochromatic Gaussian beam source.

            This source models a diffraction limited Gaussian beam characterized by its
            wavelength, polarization, numerical aperture, and optical power.

            Notes
            -----
            The numerical_aperture is treated as dimensionless and is used to derive the beam waist and
            related quantities in the C++ implementation.
        )pbdoc"
        )
        .def(
            py::init([](
                const py::object& wavelength,
                const py::object& polarization,
                const py::object& numerical_aperture,
                const py::object& optical_power
            ) {
                return std::make_shared<Gaussian>(
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>(),
                    Casting::cast_py_to_polarization_state(polarization),
                    numerical_aperture.cast<double>(),
                    optical_power.attr("to")("watt").attr("magnitude").cast<double>()
                );
            }),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("numerical_aperture"),
            py::arg("optical_power"),
            R"pbdoc(
                Construct a Gaussian beam source.

                Parameters
                ----------
                wavelength : pint.Quantity
                    Wavelength in vacuum. Must be convertible to meters.
                polarization : object
                    Polarization specification. See :class:`~PyMieSim.single.source.PlaneWave`
                    for accepted formats.
                numerical_aperture : pint.Quantity
                    Numerical aperture. Must be convertible to a dimensionless quantity.
                optical_power : pint.Quantity
                    Optical power. Must be convertible to watts.

                Returns
                -------
                Gaussian
                    Configured Gaussian source.
            )pbdoc"
        )
        .def_property_readonly(
            "waist",
            [ureg](const Gaussian& self) {
                return (py::float_(self.waist) * ureg.attr("meter")).attr("to_compact")();
            },
            R"pbdoc(
                Beam waist radius.

                Returns
                -------
                pint.Quantity
                    Waist radius with units of length.
            )pbdoc"
        )
        .def_property_readonly(
            "peak_intensity",
            [ureg](const Gaussian& self) {
                return (py::float_(self.peak_intensity) * ureg.attr("W/meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Peak on axis intensity.

                Returns
                -------
                pint.Quantity
                    Peak intensity with units of power per area.
            )pbdoc"
        )
        .def_property_readonly(
            "area",
            [ureg](const Gaussian& self) {
                return (py::float_(self.area) * ureg.attr("meter**2")).attr("to_compact")();
            },
            R"pbdoc(
                Effective beam area.

                Returns
                -------
                pint.Quantity
                    Area with units of square length.
            )pbdoc"
        )
        .def_property_readonly(
            "numerical_aperture",
            [](const Gaussian& self) {
                return py::float_(self.numerical_aperture);
            },
            R"pbdoc(
                Numerical aperture.

                Returns
                -------
                pint.Quantity
                    Dimensionless numerical aperture.
            )pbdoc"
        )
        .def_property_readonly(
            "optical_power",
            [ureg](const Gaussian& self) {
                return (py::float_(self.optical_power) * ureg.attr("watt")).attr("to_compact")();
            },
            R"pbdoc(
                Optical power.

                Returns
                -------
                pint.Quantity
                    Optical power with units of watts.
            )pbdoc"
        )
        ;

}