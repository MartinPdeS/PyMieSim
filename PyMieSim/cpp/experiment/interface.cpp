#include <pybind11/pybind11.h>


#include "experiment/setup/interface.cpp"
#include "experiment/sets/interface.cpp"

namespace py = pybind11;

PYBIND11_MODULE(interface_experiment, module) {
    module.doc() = R"pbdoc(
        Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.

        This module provides C++ bindings for the PyMieSim Python package, which implements the Lorenz-Mie Theory (LMT) for light scattering by spherical particles and other scatterers.
    )pbdoc";

    register_sets(module);

    register_setup(module);
}
