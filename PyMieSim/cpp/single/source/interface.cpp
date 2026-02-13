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

static py::array_t<std::complex<double>> first_row_to_numpy(const JonesVector& polarization)
{
    const auto& rows = polarization.elements();
    if (rows.empty()) {
        throw std::runtime_error("polarization has no elements.");
    }

    py::array_t<std::complex<double>> out(2);
    auto w = out.mutable_unchecked<1>();

    w(0) = rows[0][0];
    w(1) = rows[0][1];

    return out;
}

static py::array_t<std::complex<double>> element_to_numpy_2d(const JonesVector& polarization)
{
    const auto& rows = polarization.elements();
    const py::ssize_t n_rows = static_cast<py::ssize_t>(rows.size());

    py::array_t<std::complex<double>> out({n_rows, static_cast<py::ssize_t>(2)});
    auto w = out.mutable_unchecked<2>();

    for (py::ssize_t i = 0; i < n_rows; ++i) {
        w(i, 0) = rows[static_cast<std::size_t>(i)][0];
        w(i, 1) = rows[static_cast<std::size_t>(i)][1];
    }

    return out;
}


static JonesVector to_cpp_polarization(py::object polarization_like)
{
    // Case 1: already a bound C++ JonesVector (or derived) object
    if (py::isinstance<JonesVector>(polarization_like)) {
        return polarization_like.cast<JonesVector>();
    }

    // Case 2 and 3: try to normalize via the existing Python polarization module
    py::object polarization_mod = py::module_::import("PyMieSim.single.polarization");
    py::object BasePolarization = polarization_mod.attr("BasePolarization");
    py::object Linear = polarization_mod.attr("Linear");

    if (!py::isinstance(polarization_like, BasePolarization)) {
        polarization_like = Linear(polarization_like);
    }

    // Your Python classes store the Jones matrix in .element (shape (N,2))
    py::object element_obj = polarization_like.attr("element");

    py::array element_any = py::array::ensure(element_obj);
    if (!element_any) {
        throw std::invalid_argument("polarization.element must be array like.");
    }

    py::array_t<std::complex<double>, py::array::c_style | py::array::forcecast> element_arr(element_any);

    if (element_arr.ndim() != 2) {
        throw std::invalid_argument("polarization.element must be a 2D array with shape (N, 2).");
    }
    if (element_arr.shape(1) != 2) {
        throw std::invalid_argument("polarization.element must have second dimension == 2.");
    }
    if (element_arr.shape(0) < 1) {
        throw std::invalid_argument("polarization.element must have at least one row.");
    }

    const std::size_t n_rows = static_cast<std::size_t>(element_arr.shape(0));
    auto r = element_arr.unchecked<2>();

    JonesVector::Element rows;
    rows.reserve(n_rows);

    for (std::size_t i = 0; i < n_rows; ++i) {
        rows.push_back(JonesVector::Row{
            r(static_cast<py::ssize_t>(i), 0),
            r(static_cast<py::ssize_t>(i), 1),
        });
    }

    return JonesVector(std::move(rows));
}




void register_sources(py::module_& module)
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
                return (py::float_(self.wavelength) * ureg.attr("meter")).attr("to_compact")();
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
        .def_property_readonly(
            "polarization_element",
            [ureg](const BaseSource& self) {
                return element_to_numpy_2d(self.polarization);
            },
            R"pbdoc(
                Polarization Jones matrix.

                The polarization is returned as a complex array with shape ``(N, 2)``.
                Each row corresponds to a Jones vector ``[Ex, Ey]``.

                Returns
                -------
                numpy.ndarray
                    Complex array of shape ``(N, 2)``.
            )pbdoc"
        )
        .def_property_readonly(
            "jones_vector",
            [ureg](const BaseSource& self) {
                return first_row_to_numpy(self.polarization);
            },
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
                py::object polarization,
                py::object amplitude
            ) {
                py::object units_mod = py::module_::import("PyMieSim.units");
                py::object UnitLength = units_mod.attr("Length");
                py::object UnitAmplitude = units_mod.attr("ElectricField");

                polarization = py::reinterpret_borrow<py::object>(polarization);
                JonesVector polarization_cpp = to_cpp_polarization(polarization);

                amplitude = UnitAmplitude.attr("check")(amplitude);
                wavelength = UnitLength.attr("check")(wavelength);

                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                const double amplitude_v_per_m =
                    amplitude.attr("to")(ureg.attr("volt/meter")).attr("magnitude").cast<double>();

                return std::make_shared<Planewave>(
                    wavelength_meter,
                    std::move(polarization_cpp),
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
            The NA is treated as dimensionless and is used to derive the beam waist and
            related quantities in the C++ implementation.
        )pbdoc"
        )
        .def(
            py::init([ureg](
                py::object wavelength,
                py::object polarization,
                py::object NA,
                py::object optical_power
            ) {
                py::object units_mod = py::module_::import("PyMieSim.units");
                py::object UnitLength = units_mod.attr("Length");

                polarization = py::reinterpret_borrow<py::object>(polarization);
                JonesVector polarization_cpp = to_cpp_polarization(polarization);

                wavelength = UnitLength.attr("check")(wavelength);

                const double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                const double NA_value =
                    NA.attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();

                const double optical_power_watt =
                    optical_power.attr("to")(ureg.attr("watt")).attr("magnitude").cast<double>();

                return std::make_shared<Gaussian>(
                    wavelength_meter,
                    std::move(polarization_cpp),
                    NA_value,
                    optical_power_watt
                );
            }),
            py::arg("wavelength"),
            py::arg("polarization"),
            py::arg("NA"),
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
                NA : pint.Quantity
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
            "NA",
            [ureg](const Gaussian& self) {
                return (py::float_(self.NA) * ureg.attr("dimensionless")).attr("to_compact")();
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