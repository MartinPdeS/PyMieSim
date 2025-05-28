#include <pybind11/pybind11.h>
#include "fibonacci/fibonacci.h"

namespace py = pybind11;

PYBIND11_MODULE(interface_fibonacci, module) {
    module.doc() = "Generalized Lorenz-Mie Theory (GLMT) C++ binding module for light scattering from a spherical scatterer";

}
