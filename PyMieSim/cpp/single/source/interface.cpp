#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "source.h"
#include "pint/pint.h"

namespace py = pybind11;

py::object ureg = get_shared_ureg();

void register_sources(py::module_& module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<BaseSource, std::shared_ptr<BaseSource>>(module, "BindedBaseSource")
        .def(py::init<>());

    py::class_<Planewave, BaseSource, std::shared_ptr<Planewave>>(module, "PLANEWAVE")
        .def(
            py::init([](py::object wavelength, std::vector<complex128> jones_vector, py::object amplitude) {

                double wavelength_meter = wavelength.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                double amplitude_vpm = amplitude.attr("to")(ureg.attr("volt/meter")).attr("magnitude").cast<double>();

                return std::make_shared<Planewave>(
                    wavelength_meter,
                    std::move(jones_vector),
                    amplitude_vpm
                );
            }),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            "Constructs a Planewave source with specified optical properties."
        )
        .def_property_readonly(
            "wavelength",
            [&](const Planewave &self) {
                return (py::float_(self.wavelength) * ureg.attr("meter")).attr("to_compact")();
            },
            "Wavelength of the source."
        )
        .def_readonly(
            "_cpp_jones_vector",
            &Planewave::jones_vector,
            "Jones vector of the source."
        )
        .def_property_readonly(
            "amplitude",
            [](const Planewave &self) {
                return (py::float_(self.amplitude) * ureg.attr("volt/meter")).attr("to_compact")();
            },
            "Electric field amplitude of the source."
        )
        .def_property_readonly(
            "wavenumber",
            [](const Planewave &self) {
                return (py::float_(self.wavenumber) * ureg.attr("1/meter")).attr("to_compact")();
            },
            "Wavenumber of the source."
        )
        ;

    py::class_<Gaussian, BaseSource, std::shared_ptr<Gaussian>>(module, "GAUSSIAN")
        .def(
            py::init([](py::object wavelength, std::vector<complex128> jones_vector, py::object NA, py::object optical_power) {

                double wavelength_meter =
                    wavelength.attr("to")(ureg.attr("meter")).attr("magnitude").cast<double>();

                double NA_value =
                    (NA * ureg.attr("dimensionless")).attr("to")(ureg.attr("dimensionless")).attr("magnitude").cast<double>();

                double optical_power_watt =
                    optical_power.attr("to")(ureg.attr("watt")).attr("magnitude").cast<double>();

                return std::make_shared<Gaussian>(
                    wavelength_meter,
                    std::move(jones_vector),
                    NA_value,
                    optical_power_watt
                );
            }),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            "Constructs a Gaussian source with specified optical properties."
        )
        .def_property_readonly(
            "wavelength",
            [&](const Gaussian &self) {
                return (py::float_(self.wavelength) * ureg.attr("meter")).attr("to_compact")();
            },
            "Wavelength of the source."
        )
        .def_property_readonly(
            "wavenumber",
            [](const Gaussian &self) {
                return (py::float_(self.wavenumber) * ureg.attr("1/meter")).attr("to_compact")();
            },
            "Wavenumber of the source."
        )
        .def_property_readonly(
            "waist",
            [](const Gaussian &self) {
                return (py::float_(self.waist) * ureg.attr("meter")).attr("to_compact")();
            },
            "Beam waist of the Gaussian source."
        )
        .def_property_readonly(
            "peak_intensity",
            [](const Gaussian &self) {
                return (py::float_(self.peak_intensity) * ureg.attr("W/meter**2")).attr("to_compact")();
            },
            "Peak intensity of the Gaussian source."
        )
        .def_property_readonly(
            "area",
            [](const Gaussian &self) {
                return (py::float_(self.area) * ureg.attr("meter**2")).attr("to_compact")();
            },
            "Beam area of the Gaussian source."
        )
        .def_readonly(
            "_cpp_jones_vector",
            &Gaussian::jones_vector,
            "Jones vector of the source."
        )
        .def_property_readonly(
            "amplitude",
            [](const Gaussian &self) {
                return (py::float_(self.amplitude) * ureg.attr("volt/meter")).attr("to_compact")();
            },
            "Electric field amplitude of the source."
        )
        .def_property_readonly(
            "NA",
            [](const Gaussian &self) {
                return (py::float_(self.NA) * ureg.attr("dimensionless")).attr("to_compact")();
            },
            "Numerical Aperture of the source."
        )
        .def_property_readonly(
            "optical_power",
            [](const Gaussian &self) {
                return (py::float_(self.optical_power) * ureg.attr("watt")).attr("to_compact")();
            },
            "Optical power of the source."
        );
}
