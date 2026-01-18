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

#include "source.h"
#include "pint/pint.h"

namespace py = pybind11;

py::object ureg = get_shared_ureg();

// Accept polarization as:
// 1) a bound C++ polarization object (JonesVector / Linear / RightCircular / LeftCircular) from our module
// 2) a legacy Python polarization object with .element shaped (N,2)
// 3) anything else: wrap with Python Linear(...) if possible (for backward compatibility)
//
// Ultimately returns a C++ JonesVector with element rows.
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

static py::array_t<std::complex<double>> first_row_to_numpy(const JonesVector& polarization)
{
    const auto& rows = polarization.elements();
    if (rows.empty()) {
        throw std::runtime_error("polarization has no elements.");
    }

    py::array_t<std::complex<double>> out({static_cast<py::ssize_t>(2)});
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


void register_sources(py::module_& module)
{
    module.doc() = "Source bindings for PyMieSim.";

    py::class_<BaseSource, std::shared_ptr<BaseSource>>(module, "BaseSource")
        .def(py::init<>())
        .def_property_readonly(
            "wavelength",
            [](const BaseSource& self) {
                return (py::float_(self.wavelength) * ureg.attr("meter")).attr("to_compact")();
            }
        )
        .def_property_readonly(
            "wavenumber",
            [](const BaseSource& self) {
                return (py::float_(self.wavenumber) * ureg.attr("1/meter")).attr("to_compact")();
            }
        )
        .def_property_readonly(
            "amplitude",
            [](const BaseSource& self) {
                return (py::float_(self.amplitude) * ureg.attr("volt/meter")).attr("to_compact")();
            }
        )
        // Expose polarization as a 2D complex array (N,2) plus a convenience first row.
        .def_property_readonly(
            "polarization_element",
            [](const BaseSource& self) {
                return element_to_numpy_2d(self.polarization);
            },
            "Polarization Jones matrix as a complex array of shape (N, 2)."
        )
        .def_property_readonly(
            "jones_vector",
            [](const BaseSource& self) {
                return first_row_to_numpy(self.polarization);
            },
            "First Jones vector row [Ex, Ey] as a complex array of shape (2,)."
        );

    py::class_<Planewave, BaseSource, std::shared_ptr<Planewave>>(module, "PlaneWave")
        .def(
            py::init([](
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
            "Constructs a PlaneWave source with specified optical properties."
        );

    py::class_<Gaussian, BaseSource, std::shared_ptr<Gaussian>>(module, "Gaussian")
        .def(
            py::init([](
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
            "Constructs a Gaussian source with specified optical properties."
        )
        .def_property_readonly(
            "waist",
            [](const Gaussian& self) {
                return (py::float_(self.waist) * ureg.attr("meter")).attr("to_compact")();
            }
        )
        .def_property_readonly(
            "peak_intensity",
            [](const Gaussian& self) {
                return (py::float_(self.peak_intensity) * ureg.attr("W/meter**2")).attr("to_compact")();
            }
        )
        .def_property_readonly(
            "area",
            [](const Gaussian& self) {
                return (py::float_(self.area) * ureg.attr("meter**2")).attr("to_compact")();
            }
        )
        .def_property_readonly(
            "NA",
            [](const Gaussian& self) {
                return (py::float_(self.NA) * ureg.attr("dimensionless")).attr("to_compact")();
            }
        )
        .def_property_readonly(
            "optical_power",
            [](const Gaussian& self) {
                return (py::float_(self.optical_power) * ureg.attr("watt")).attr("to_compact")();
            }
        )
        ;
    }