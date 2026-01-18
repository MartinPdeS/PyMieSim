#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(interface_pint, module) {
    module.def(
        "set_ureg",
        [module](py::object ureg_object) mutable {
            py::gil_scoped_acquire gil;
            if (ureg_object.is_none()) {
                throw std::runtime_error("set_ureg: ureg is None.");
            }
            module.attr("_ureg") = std::move(ureg_object);
        },
        py::arg("ureg_object")
    );

    module.def(
        "get_ureg",
        [module]() mutable {
            py::gil_scoped_acquire gil;
            if (!py::hasattr(module, "_ureg")) {
                throw std::runtime_error("get_ureg: ureg not initialized. Call set_ureg first.");
            }
            py::object ureg = module.attr("_ureg");
            if (ureg.is_none()) {
                throw std::runtime_error("get_ureg: ureg is None. Call set_ureg with a valid registry.");
            }
            return ureg;
        }
    );
}
