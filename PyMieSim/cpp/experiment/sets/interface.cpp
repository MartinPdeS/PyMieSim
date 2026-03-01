#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers
#include <pybind11/complex.h> // For std::complex support
#include <limits>

#include "properties_set/interface.h"
#include "source_set/interface.h"
#include "scatterer_set/interface.h"
#include "detector_set/interface.h"

namespace py = pybind11;

void register_sets(py::module& module) {
    register_properties_set(module);
    register_source_set(module);
    register_scatterer_set(module);
    register_detector_set(module);
}

// -
