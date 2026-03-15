#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "./setup.cpp"
#include <utils/numpy_interface.h>
#include <utils/defines.h>

namespace py = pybind11;

#define DEFINE_GETTER_INTERFACE(property) \
    .def("get_"  #property, \
        [](Setup& self) { \
            auto [output, shape] = self.get_data<&BaseScatterer::get_##property>(); \
            return vector_move_from_numpy(output, shape); \
        }, \
        "Retrieves the scattering property for a scatterer" \
    ) \
    .def("get_" #property "_sequential", \
        [](Setup& self) { \
            std::vector<double> output = self.get_data_sequential<&BaseScatterer::get_##property>(); \
            return vector_move_from_numpy(output, {output.size()}); \
        }, \
        "Retrieves the scattering property for a scatterer" \
    )

#define DEFINE_GETTERS_INTERFACE_1(A) \
    DEFINE_GETTER_INTERFACE(A)

#define DEFINE_GETTERS_INTERFACE_2(A, B) \
    DEFINE_GETTER_INTERFACE(A) \
    DEFINE_GETTER_INTERFACE(B)

#define DEFINE_GETTERS_INTERFACE_3(A, B, C) \
    DEFINE_GETTER_INTERFACE(A) \
    DEFINE_GETTER_INTERFACE(B) \
    DEFINE_GETTER_INTERFACE(C)

#define DEFINE_GETTERS_INTERFACE_4(A, B, C, D) \
    DEFINE_GETTER_INTERFACE(A) \
    DEFINE_GETTER_INTERFACE(B) \
    DEFINE_GETTER_INTERFACE(C) \
    DEFINE_GETTER_INTERFACE(D)

#define DEFINE_GETTERS_INTERFACE_6(A, B, C, D, E, F) \
    DEFINE_GETTER_INTERFACE(A) \
    DEFINE_GETTER_INTERFACE(B) \
    DEFINE_GETTER_INTERFACE(C) \
    DEFINE_GETTER_INTERFACE(D) \
    DEFINE_GETTER_INTERFACE(E) \
    DEFINE_GETTER_INTERFACE(F)

#define DEFINE_GETTERS_INTERFACE_7(A, B, C, D, E, F, G) \
    DEFINE_GETTER_INTERFACE(A) \
    DEFINE_GETTER_INTERFACE(B) \
    DEFINE_GETTER_INTERFACE(C) \
    DEFINE_GETTER_INTERFACE(D) \
    DEFINE_GETTER_INTERFACE(E) \
    DEFINE_GETTER_INTERFACE(F) \
    DEFINE_GETTER_INTERFACE(G)


PYBIND11_MODULE(_setup, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    pybind11::class_<Setup>(module, "Setup")
        .def(
            py::init<
                std::shared_ptr<ScattererSet>,
                std::shared_ptr<BaseSourceSet>,
                std::shared_ptr<BaseDetectorSet>,
                bool
            >(),
            py::arg("scatterer_set"),
            py::arg("source_set"),
            py::arg("detector_set") = nullptr,
            py::arg("debug_mode") = false,
            R"pbdoc(
                Setup class for conducting Lorenz-Mie Theory (LMT) simulations.

                This class manages the interaction between scatterers, sources, and detectors,
                providing methods to compute scattering properties, coupling coefficients,
                and far-field patterns.

                Parameters
                ----------
                debug_mode : bool, optional
                    If set to True, enables debug printing for tracing computations. Default is True.
            )pbdoc"
        )
        .def_readwrite("scatterer_set", &Setup::scatterer_set)
        .def_readwrite("source_set", &Setup::source_set)
        .def_readwrite("detector_set", &Setup::detector_set)
        .def_readonly("debug_mode", &Setup::debug_mode)
        .def_readonly("array_shape", &Setup::array_shape)
        .def_readonly("total_iterations", &Setup::total_iterations)
        .def(
            "get_coupling_sequential",
            [](Setup& self) {

                std::vector<double> coupling_array = self.get_coupling_sequential();
                return vector_move_from_numpy(coupling_array, {coupling_array.size()});
            },
            R"pbdoc(
                Retrieves the coupling power for a combination of scatterers, sources, and detectors in a sequential manner.
                Returns
                -------
                numpy.ndarray
                    A numpy array containing the coupling coefficients.
            )pbdoc"
        )
        .def(
            "get_coupling",
            [](Setup& self) {

                auto [coupling_array, coupling_shape] = self.get_coupling();
                return vector_move_from_numpy(coupling_array, coupling_shape);
            },
            R"pbdoc(
                Retrieves the coupling power for a combination of scatterers, sources, and detectors.
            )pbdoc"
        )
        .def(
            "_get_farfields",
            [](Setup& self, const FibonacciMesh& mesh, const double distance){

                auto [farfield_array, farfield_shape] = self.get_farfields(mesh, distance);
                return vector_move_from_numpy(farfield_array, farfield_shape);
            },
            pybind11::arg("mesh"),
            pybind11::arg("distance") = 1.0,
            R"pbdoc(
                Retrieves the far-field patterns for a combination of scatterers and sources on a Fibonacci mesh.

                Parameters
                ----------
                mesh : FibonacciMesh
                    The Fibonacci mesh for far-field sampling.
                distance : float, optional
                    The distance at which to compute the far-fields. Default is 1.0.

                Returns
                -------
                numpy.ndarray
                    A numpy array containing the far-field patterns.
            )pbdoc"
        )

        DEFINE_GETTERS_INTERFACE_6(a1, a2, a3, b1, b2, b3)

        DEFINE_GETTERS_INTERFACE_6(a11, a12, a13, b11, b12, b13)
        DEFINE_GETTERS_INTERFACE_6(a21, a22, a23, b21, b22, b23)

        DEFINE_GETTERS_INTERFACE_7(Qsca, Qext, Qabs, Qpr, Qforward, Qback, Qratio)
        DEFINE_GETTERS_INTERFACE_7(Csca, Cext, Cabs, Cpr, Cforward, Cback, Cratio)
        DEFINE_GETTERS_INTERFACE_1(g)
        ;
}
