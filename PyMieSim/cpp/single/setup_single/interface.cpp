#include <pybind11/pybind11.h>
#include "./setup_single.h"
#include <utils/numpy_interface.h>
#include <utils/defines.h>

namespace py = pybind11;



void register_setup_single(py::module& module) {
    py::class_<Setup, std::shared_ptr<Setup>>(module,
        "SETUP_SINGLE",
        py::module_local(),
            R"pbdoc(
                Setup for single scatterer Lorenz-Mie simulations.

                This class encapsulates the scatterer, source, and detector
                configurations for a single scattering simulation.
            )pbdoc"
        )
        .def(
            py::init(
                [](
                    std::shared_ptr<BaseScatterer> scatterer,
                    std::shared_ptr<BaseSource> source,
                    std::shared_ptr<Photodiode> detector,
                    bool debug_mode = false
                ) {
                    return std::make_shared<Setup>(scatterer, source, detector, debug_mode);
                }
            ),
            py::arg("scatterer"),
            py::arg("source"),
            py::arg("detector"),
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
        ;


}
