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
    )pbdoc")
        .def(py::init<>(),
        R"pbdoc(
            Construct a source base instance.

            Notes
            -----
            This constructor is primarily used for binding completeness.
            In user code, you typically instantiate a concrete source such as
            :class:`~PyMieSim.single.source.PlaneWave` or
            :class:`~PyMieSim.single.source.Gaussian`.
        )pbdoc")
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
        );


    py::class_<Planewave, BaseSource, std::shared_ptr<Planewave>>(module, "PlaneWave",
        R"pbdoc(
            Monochromatic plane wave source.

            This source models a spatially uniform incident field with a specified
            wavelength, polarization, and electric field amplitude.
        )pbdoc")
        .def(
            py::init([ureg](
                py::object wavelength,
                PolarizationState polarization,
                py::object amplitude
            ) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                const double amplitude_v_per_m =
                    amplitude.attr("to")("volt/meter").attr("magnitude").cast<double>();

                return std::make_shared<Planewave>(
                    wavelength_meter,
                    polarization,
                    amplitude_v_per_m
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
            py::init([ureg](
                py::object wavelength,
                PolarizationState polarization,
                py::object numerical_aperture,
                py::object optical_power
            ) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                const double numerical_aperture_value =
                    numerical_aperture.attr("to")("dimensionless").attr("magnitude").cast<double>();

                const double optical_power_watt =
                    optical_power.attr("to")("watt").attr("magnitude").cast<double>();

                return std::make_shared<Gaussian>(
                    wavelength_meter,
                    polarization,
                    numerical_aperture_value,
                    optical_power_watt
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
            [ureg](const Gaussian& self) {
                return (py::float_(self.numerical_aperture) * ureg.attr("dimensionless")).attr("to_compact")();
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
        );

}