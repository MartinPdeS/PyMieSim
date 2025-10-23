#include <pybind11/pybind11.h>
#include "experiment/experiment.cpp"

#define EFFICIENCY_PROPERTY(scatterer, name) \
    .def("get_" #scatterer "_Q" #name, &Experiment::get_##scatterer##_Q##name, py::arg("source_set"), py::arg("detector_set"), "Retrieves the scattering efficiency (Q" #name ") for a" #scatterer) \
    .def("get_" #scatterer "_C" #name, &Experiment::get_##scatterer##_C##name, py::arg("source_set"), py::arg("detector_set"), "Retrieves the scattering cross-section (C" #name ") for a" #scatterer) \
    .def("get_" #scatterer "_Q" #name "_sequential", &Experiment::get_##scatterer##_Q##name##_sequential, py::arg("source_set"), py::arg("detector_set"), "Retrieves the scattering efficiency (Q" #name ") for a" #scatterer " in a sequential manner.") \
    .def("get_" #scatterer "_C" #name "_sequential", &Experiment::get_##scatterer##_C##name##_sequential, py::arg("source_set"), py::arg("detector_set"), "Retrieves the scattering cross-section (C" #name ") for a" #scatterer " in a sequential manner.")

#define COEFFICIENT_PROPERTY(scatterer, number) \
    .def("get_" #scatterer "_a" #number, &Experiment::get_##scatterer##_a##number, py::arg("source_set"), py::arg("detector_set"), "Retrieves the electric multipole coefficient (a" #number ") coefficient for a " #scatterer) \
    .def("get_" #scatterer "_b" #number, &Experiment::get_##scatterer##_b##number, py::arg("source_set"), py::arg("detector_set"), "Retrieves the magnetic multipole coefficient (b" #number ") coefficient for a " #scatterer) \
    .def("get_" #scatterer "_a" #number "_sequential", &Experiment::get_##scatterer##_a##number##_sequential, py::arg("source_set"), py::arg("detector_set"), "Retrieves the electric multipole coefficient (a" #number ") coefficient for a " #scatterer " in a sequential manner.") \
    .def("get_" #scatterer "_b" #number "_sequential", &Experiment::get_##scatterer##_b##number##_sequential, py::arg("source_set"), py::arg("detector_set"), "Retrieves the magnetic multipole coefficient (b" #number ") coefficient for a " #scatterer " in a sequential manner.")

#define COUPLING_AND_G(scatterer) \
     .def("get_" #scatterer "_g", &Experiment::get_##scatterer##_g, py::arg("source_set"), py::arg("detector_set"), "Retrieves the asymmetry parameter (g) for a " #scatterer) \
     .def("get_" #scatterer "_g_sequential", &Experiment::get_##scatterer##_g_sequential, py::arg("source_set"), py::arg("detector_set"), "Retrieves the asymmetry parameter (g) for a " #scatterer " in a sequential manner.")



#define DEFINE_GETTERS(property) \
    .def( \
        "get_"  #property, \
        &Experiment::get_data<&BaseScatterer::get_##property>, \
        py::arg("scatterer_set"), \
        py::arg("source_set"), \
        py::arg("detector_set"), \
        "Retrieves the scattering property for a scatterer") \
    .def( \
        "get_" #property "_sequential", \
        &Experiment::get_data_sequential<&BaseScatterer::get_##property>, \
        py::arg("scatterer_set"), \
        py::arg("source_set"), \
        py::arg("detector_set"), \
        "Retrieves the scattering property for a scatterer")

namespace py = pybind11;

PYBIND11_MODULE(interface_experiment, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    py::class_<Experiment>(module, "EXPERIMENT")
        .def(py::init<bool>(), py::arg("debug_mode") = true, "Constructs an Experiment object.")

        .def("get_coupling_sequential", &Experiment::get_coupling_sequential, py::arg("scatterer_set"), py::arg("source_set"), py::arg("detector_set"), "Retrieves the coupling power for a scatterer")
        .def("get_coupling", &Experiment::get_coupling, py::arg("scatterer_set"), py::arg("source_set"), py::arg("detector_set"), "Retrieves the coupling power for a scatterer")
        .def("_get_farfields", &Experiment::get_farfields, py::arg("scatterer_set"), py::arg("source_set"), py::arg("mesh"), py::arg("distance") = 1, "Retrieves the farfields for a scatterer")

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
