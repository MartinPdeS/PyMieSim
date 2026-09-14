#include "distributions.h"

#include <pybind11/stl.h>

PYBIND11_MODULE(distributions, module) {
    module.doc() = "Native particle-size distributions for independent-scattering averages.";

    module.def("legendre", &legendre, py::arg("order"), "Return Gauss-Legendre nodes and weights on [-1, 1].");

    module.def("hermite", &hermite, py::arg("order"), "Return Gauss-Hermite nodes and weights for exp(-x**2).");

    py::class_<ParticleSizeDistribution>(module, "ParticleSizeDistribution",
                                         R"pbdoc(
            Discrete particle diameters and normalized number fractions.

            Optical averaging uses the non-interacting, independent-scattering
            approximation. Refer to PackLab when particle correlations matter.
        )pbdoc")
        .def(py::init<py::object, py::object>(), py::arg("diameters"), py::arg("number_weights"),
             "Create a discrete distribution from unit-aware diameters and number weights.")
        .def_property_readonly("diameters", &ParticleSizeDistribution::diameters,
                               "Diameter nodes as a Pint quantity in meters.")
        .def_property_readonly("number_fractions", &ParticleSizeDistribution::number_fractions,
                               "Normalized particle-number fractions as a NumPy array.")
        .def_static("monodisperse", &ParticleSizeDistribution::monodisperse, py::arg("diameter"),
                    "Create a distribution containing one diameter.")
        .def_static(
            "uniform",
            [](py::object minimum, py::object maximum, py::object sampling) {
                return ParticleSizeDistribution::uniform(std::move(minimum), std::move(maximum),
                                                         parse_sampling(std::move(sampling)));
            },
            py::arg("minimum_diameter"), py::arg("maximum_diameter"), py::arg("sampling") = 32)
        .def_static(
            "lognormal",
            [](py::object median, py::object geometric_std, py::object sampling) {
                py::module_ numpy = py::module_::import("numpy");
                py::object value = numpy.attr("asarray")(geometric_std, numpy.attr("float64"));
                if (value.attr("ndim").cast<int>() != 0) {
                    throw py::value_error("geometric_std must be finite and at least one");
                }
                return ParticleSizeDistribution::lognormal(std::move(median), value.attr("item")().cast<double>(),
                                                           parse_sampling(std::move(sampling)));
            },
            py::arg("median_diameter"), py::arg("geometric_std"), py::arg("sampling") = 32)
        .def_static(
            "truncated_normal",
            [](py::object mean, py::object width, py::object minimum, py::object maximum, py::object sampling) {
                return ParticleSizeDistribution::truncated_normal(std::move(mean), std::move(width), std::move(minimum),
                                                                  std::move(maximum),
                                                                  parse_sampling(std::move(sampling)));
            },
            py::arg("mean_diameter"), py::arg("standard_deviation"), py::kw_only(), py::arg("minimum_diameter"),
            py::arg("maximum_diameter"), py::arg("sampling") = 64)
        .def_static(
            "triangular",
            [](py::object minimum, py::object mode, py::object maximum, py::object sampling) {
                return ParticleSizeDistribution::triangular(std::move(minimum), std::move(mode), std::move(maximum),
                                                            parse_sampling(std::move(sampling)));
            },
            py::arg("minimum_diameter"), py::arg("mode_diameter"), py::arg("maximum_diameter"),
            py::arg("sampling") = 16)
        .def_static("mixture", &ParticleSizeDistribution::mixture, py::arg("distributions"), py::arg("number_weights"),
                    "Combine distributions using component number fractions.")
        .def("__repr__", &ParticleSizeDistribution::repr);
}
