#include <pybind11/pybind11.h>
#include "./base_scatterer.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>

namespace py = pybind11;

void register_base_scatterer(py::module_& module) {

    py::object ureg = get_shared_ureg();

    py::class_<BaseScatterer, std::shared_ptr<BaseScatterer>>(module, "BaseScatterer")
        .def(
            "init",
            &BaseScatterer::init,
            py::arg("source"),
            py::arg("max_order") = 0,
            R"pbdoc(
                Initializes the scatterer by computing size parameter, cross section, and scattering coefficients.

                Parameters
                ----------
                source : BaseSource
                    The source of the electromagnetic field, used to determine incident field properties for scatterer initialization.
                max_order : int, optional
                    The maximum order of the coefficients to compute (default is 0, which means it will be computed based on the size parameter).
            )pbdoc"
        )
        .def(
            "compute_cn_dn",
            &BaseScatterer::compute_cn_dn,
            py::arg("max_order") = 0,
            R"pbdoc(
                Computes the coefficients cn and dn for the scatterer.

                Parameters
                ----------
                max_order : int, optional
                    The maximum order of the coefficients to compute (default is 0, which means it will be computed based on the size parameter).
            )pbdoc"
        )
        .def(
            "compute_an_bn",
            &BaseScatterer::compute_an_bn,
            py::arg("max_order") = 0,
            R"pbdoc(
                Computes the coefficients an and bn for the scatterer.
                Parameters
                ----------
                max_order : int, optional
                    The maximum order of the coefficients to compute (default is 0, which means it will be computed based on the size parameter).
            )pbdoc"
        )
        .def(
            "get_coefficient",
            &BaseScatterer::get_coefficient,
            py::arg("type"),
            py::arg("order"),
            R"pbdoc(
                Returns a specific scattering coefficient for a given type and order.

                Parameters
                ----------
                type : str
                    The type of coefficient ('an', 'bn', 'cn', 'dn').
                order : int
                    The order of the scattering coefficient.

                Returns
                -------
                float
                    The value of the specified scattering coefficient.
            )pbdoc"
        )
        ;
}
