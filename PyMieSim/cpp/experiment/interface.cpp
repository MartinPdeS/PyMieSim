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

        // Setup methods
        .def("set_Sphere", &Experiment::set_sphere, "Defines a spherical scatterer for the experiment.")
        .def("set_Cylinder", &Experiment::set_cylinder, "Defines a cylindrical scatterer for the experiment.")
        .def("set_CoreShell", &Experiment::set_coreshell, "Defines a core-shell scatterer for the experiment.")

        .def("get_coupling_sequential", &Experiment::get_coupling_sequential, py::arg("scatterer_set"), py::arg("source_set"), py::arg("detector_set"), "Retrieves the coupling power for a scatterer")
        .def("get_coupling", &Experiment::get_coupling, py::arg("scatterer_set"), py::arg("source_set"), py::arg("detector_set"), "Retrieves the coupling power for a scatterer")


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



        // Sphere metrics
        EFFICIENCY_PROPERTY(Sphere, sca)  // Efficiencies
        EFFICIENCY_PROPERTY(Sphere, ext)
        EFFICIENCY_PROPERTY(Sphere, abs)
        EFFICIENCY_PROPERTY(Sphere, pr)
        EFFICIENCY_PROPERTY(Sphere, forward)
        EFFICIENCY_PROPERTY(Sphere, back)
        EFFICIENCY_PROPERTY(Sphere, ratio)
        COEFFICIENT_PROPERTY(Sphere, 1)   // Coefficients
        COEFFICIENT_PROPERTY(Sphere, 2)
        COEFFICIENT_PROPERTY(Sphere, 3)


        // Cylinder metrics
        EFFICIENCY_PROPERTY(Cylinder, sca)  // Efficiencies
        EFFICIENCY_PROPERTY(Cylinder, ext)
        EFFICIENCY_PROPERTY(Cylinder, abs)
        COEFFICIENT_PROPERTY(Cylinder, 11)  // Coefficients
        COEFFICIENT_PROPERTY(Cylinder, 21)
        COEFFICIENT_PROPERTY(Cylinder, 12)
        COEFFICIENT_PROPERTY(Cylinder, 22)
        COEFFICIENT_PROPERTY(Cylinder, 13)
        COEFFICIENT_PROPERTY(Cylinder, 23)


        // Coreshell metrics
        EFFICIENCY_PROPERTY(CoreShell, sca)  // Efficiencies
        EFFICIENCY_PROPERTY(CoreShell, ext)
        EFFICIENCY_PROPERTY(CoreShell, abs)
        EFFICIENCY_PROPERTY(CoreShell, pr)
        EFFICIENCY_PROPERTY(CoreShell, forward)
        EFFICIENCY_PROPERTY(CoreShell, back)
        EFFICIENCY_PROPERTY(CoreShell, ratio)
        COEFFICIENT_PROPERTY(CoreShell, 1)  // Coefficients
        COEFFICIENT_PROPERTY(CoreShell, 2)
        COEFFICIENT_PROPERTY(CoreShell, 3)
        ;
}

