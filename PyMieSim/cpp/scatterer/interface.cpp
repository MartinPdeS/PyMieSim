#include <pybind11/pybind11.h>

#include "scatterer/base_scatterer/interface.cpp"
#include "scatterer/sphere/interface.cpp"
#include "scatterer/coreshell/interface.cpp"
#include "scatterer/cylinder/interface.cpp"



namespace py = pybind11;

PYBIND11_MODULE(interface_scatterer, module) {
    module.doc() = R"pbdoc(
        Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.

        This module provides C++ bindings for the PyMieSim Python package, which implements the Lorenz-Mie Theory (LMT) for light scattering by spherical particles and other scatterers.
    )pbdoc";

    register_base_scatterer(module);

    register_sphere(module);

    register_coreshell(module);

    register_cylinder(module);


}
