#include <pybind11/pybind11.h>
#include "experiment/experiment.cpp"

#define EFFICIENCY_PROPERTY(scatterer, name) \
    .def("get_" #scatterer "_Q" #name, &Experiment::get_##scatterer##_Q##name, "Retrieves the scattering efficiency (Q" #name ") for a" #scatterer) \
    .def("get_" #scatterer "_C" #name, &Experiment::get_##scatterer##_C##name, "Retrieves the scattering cross-section (C" #name ") for a" #scatterer) \
    .def("get_" #scatterer "_Q" #name "_sequential", &Experiment::get_##scatterer##_Q##name##_sequential, "Retrieves the scattering efficiency (Q" #name ") for a" #scatterer " in a sequential manner.") \
    .def("get_" #scatterer "_C" #name "_sequential", &Experiment::get_##scatterer##_C##name##_sequential, "Retrieves the scattering cross-section (C" #name ") for a" #scatterer " in a sequential manner.")

#define COEFFICIENT_PROPERTY(scatterer, number) \
    .def("get_" #scatterer "_a" #number, &Experiment::get_##scatterer##_a##number, "Retrieves the electric multipole coefficient (a" #number ") coefficient for a " #scatterer) \
    .def("get_" #scatterer "_b" #number, &Experiment::get_##scatterer##_b##number, "Retrieves the magnetic multipole coefficient (b" #number ") coefficient for a " #scatterer) \
    .def("get_" #scatterer "_a" #number "_sequential", &Experiment::get_##scatterer##_a##number##_sequential, "Retrieves the electric multipole coefficient (a" #number ") coefficient for a " #scatterer " in a sequential manner.") \
    .def("get_" #scatterer "_b" #number "_sequential", &Experiment::get_##scatterer##_b##number##_sequential, "Retrieves the magnetic multipole coefficient (b" #number ") coefficient for a " #scatterer " in a sequential manner.")

#define COUPLING_AND_G(scatterer) \
     .def("get_" #scatterer "_g", &Experiment::get_##scatterer##_g, "Retrieves the asymmetry parameter (g) for a " #scatterer) \
     .def("get_" #scatterer "_coupling", &Experiment::get_##scatterer##_coupling, "Retrieves the coupling efficiency for a " #scatterer) \
     .def("get_" #scatterer "_g_sequential", &Experiment::get_##scatterer##_g_sequential, "Retrieves the asymmetry parameter (g) for a " #scatterer " in a sequential manner.") \
     .def("get_" #scatterer "_coupling_sequential", &Experiment::get_##scatterer##_coupling_sequential, "Retrieves the coupling efficiency for a " #scatterer " in a sequential manner.")


namespace py = pybind11;

PYBIND11_MODULE(interface_experiment, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    py::class_<Experiment>(module, "EXPERIMENT")
        .def(py::init<bool>(), py::arg("debug_mode") = true, "Constructs an Experiment object.")

        // Setup methods
        .def("set_Detector", &Experiment::set_detector, "Configures the detector for the experiment.")
        .def("set_Source", &Experiment::set_source, "Sets the light source for the experiment.")
        .def("set_Sphere", &Experiment::set_sphere, "Defines a spherical scatterer for the experiment.")
        .def("set_Cylinder", &Experiment::set_cylinder, "Defines a cylindrical scatterer for the experiment.")
        .def("set_CoreShell", &Experiment::set_coreshell, "Defines a core-shell scatterer for the experiment.")


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
        COUPLING_AND_G(Sphere)            // Coupling and g


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
        COUPLING_AND_G(Cylinder)            // Coupling and g


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
        COUPLING_AND_G(CoreShell)           // Coupling and g
        ;
}

