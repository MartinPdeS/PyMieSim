#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // For binding std::vector and similar STL containers


#include "properties_set/interface.h"
#include "polarization_set/interface.h"
#include "source_set/interface.h"
#include "scatterer_set/interface.h"
#include "detector_set/interface.h"

namespace py = pybind11;

void register_sets(py::module& module) {
    register_polarization_set(module);
    register_properties_set(module);
    register_source_set(module);
    register_scatterer_set(module);
    register_detector_set(module);
}

// -
