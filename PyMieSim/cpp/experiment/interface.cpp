#include <pybind11/pybind11.h>
#include "experiment/experiment.cpp"
#include "utils/numpy_interface.h"


// &Experiment::get_data<&BaseScatterer::get_##property>,
#define DEFINE_GETTERS(property) \
    .def("get_"  #property, \
        [](Experiment& self, const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) { \
            auto [output, shape] = self.get_data<&BaseScatterer::get_##property>(scatterer_set, source_set, detector_set); \
            return vector_move_from_numpy(output, shape); \
        }, \
        pybind11::arg("scatterer_set"), \
        pybind11::arg("source_set"), \
        pybind11::arg("detector_set"), \
        "Retrieves the scattering property for a scatterer" \
    ) \
    .def("get_" #property "_sequential", \
        [](Experiment& self, const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) { \
            std::vector<double> output = self.get_data_sequential<&BaseScatterer::get_##property>(scatterer_set, source_set, detector_set); \
            return vector_move_from_numpy(output, {output.size()}); \
        }, \
        pybind11::arg("scatterer_set"), \
        pybind11::arg("source_set"), \
        pybind11::arg("detector_set"), \
        "Retrieves the scattering property for a scatterer" \
    )

PYBIND11_MODULE(interface_experiment, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    pybind11::class_<Experiment>(module, "EXPERIMENT")
        .def(
            pybind11::init<bool>(),
            pybind11::arg("debug_mode") = true,
            R"pbdoc(
                Experiment class for conducting Lorenz-Mie Theory (LMT) simulations.

                This class manages the interaction between scatterers, sources, and detectors,
                providing methods to compute scattering properties, coupling coefficients,
                and far-field patterns.

                Parameters
                ----------
                debug_mode : bool, optional
                    If set to True, enables debug printing for tracing computations. Default is True.
            )pbdoc"
        )
        .def("get_coupling_sequential",
            [](Experiment& self, const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) {

                std::vector<double> coupling_array = self.get_coupling_sequential(scatterer_set, source_set, detector_set);
                return vector_move_from_numpy(coupling_array, {coupling_array.size()});
            },
            pybind11::arg("scatterer_set"),
            pybind11::arg("source_set"),
            pybind11::arg("detector_set"),
            R"pbdoc(
                Retrieves the coupling power for a combination of scatterers, sources, and detectors in a sequential manner.

                Parameters
                ----------
                scatterer_set : ScattererSet
                    The set of scatterers.
                source_set : BaseSourceSet
                    The set of sources.
                detector_set : DetectorSet
                    The set of detectors.

                Returns
                -------
                numpy.ndarray
                    A numpy array containing the coupling coefficients.
            )pbdoc"
        )
        .def("get_coupling",
            [](Experiment& self, const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) {

                auto [coupling_array, coupling_shape] = self.get_coupling(scatterer_set, source_set, detector_set);
                return vector_move_from_numpy(coupling_array, coupling_shape);
            },
            pybind11::arg("scatterer_set"),
            pybind11::arg("source_set"),
            pybind11::arg("detector_set"),
            R"pbdoc(
                Retrieves the coupling power for a combination of scatterers, sources, and detectors.

                Parameters
                ----------
                scatterer_set : ScattererSet
                    The set of scatterers.
                source_set : BaseSourceSet
                    The set of sources.
                detector_set : DetectorSet
                    The set of detectors.
            )pbdoc"
        )
        .def("_get_farfields",
            [](Experiment& self, const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const FibonacciMesh& mesh, const double distance){

                auto [farfield_array, farfield_shape] = self.get_farfields(scatterer_set, source_set, mesh, distance);
                return vector_move_from_numpy(farfield_array, farfield_shape);
            },
            pybind11::arg("scatterer_set"),
            pybind11::arg("source_set"),
            pybind11::arg("mesh"),
            pybind11::arg("distance") = 1,
            R"pbdoc(
                Retrieves the far-field patterns for a combination of scatterers and sources on a Fibonacci mesh.

                Parameters
                ----------
                scatterer_set : ScattererSet
                    The set of scatterers.
                source_set : BaseSourceSet
                    The set of sources.
                mesh : FibonacciMesh
                    The Fibonacci mesh for far-field sampling.
                distance : float, optional
                    The distance at which to compute the far-fields. Default is 1.

                Returns
                -------
                numpy.ndarray
                    A numpy array containing the far-field patterns.
            )pbdoc"
        )

        DEFINE_GETTERS(a1)
        DEFINE_GETTERS(a2)
        DEFINE_GETTERS(a3)
        DEFINE_GETTERS(b1)
        DEFINE_GETTERS(b2)
        DEFINE_GETTERS(b3)

        DEFINE_GETTERS(a11)
        DEFINE_GETTERS(a12)
        DEFINE_GETTERS(a13)
        DEFINE_GETTERS(b11)
        DEFINE_GETTERS(b12)
        DEFINE_GETTERS(b13)
        DEFINE_GETTERS(a21)
        DEFINE_GETTERS(a22)
        DEFINE_GETTERS(a23)
        DEFINE_GETTERS(b21)
        DEFINE_GETTERS(b22)
        DEFINE_GETTERS(b23)


        DEFINE_GETTERS(Qsca)
        DEFINE_GETTERS(Qext)
        DEFINE_GETTERS(Qext)
        DEFINE_GETTERS(Qabs)
        DEFINE_GETTERS(Qpr)
        DEFINE_GETTERS(Qforward)
        DEFINE_GETTERS(Qback)
        DEFINE_GETTERS(Qratio)

        DEFINE_GETTERS(Csca)
        DEFINE_GETTERS(Cext)
        DEFINE_GETTERS(Cext)
        DEFINE_GETTERS(Cabs)
        DEFINE_GETTERS(Cpr)
        DEFINE_GETTERS(Cforward)
        DEFINE_GETTERS(Cback)
        DEFINE_GETTERS(Cratio)

        DEFINE_GETTERS(g)
        ;
}
