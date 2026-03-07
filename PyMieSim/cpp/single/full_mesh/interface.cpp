#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "full_mesh.h"
#include <pint/pint.h>

namespace py = pybind11;

PYBIND11_MODULE(full_mesh, module)
{
    py::object ureg = get_shared_ureg();

    pybind11::class_<FullSteradian, std::shared_ptr<FullSteradian>>(module, "FullMesh")

    ;
}
