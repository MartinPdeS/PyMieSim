#include <pybind11/pybind11.h>
#include "experiment/includes/experiment.cpp"

#define EFFICIENCY_PROPERTY(scatterer, name) \
    .def("get_" #scatterer "_Q" #name, &Experiment::get_##scatterer##_Q##name, "Retrieves the scattering efficiency (Q" #name ") for a" #scatterer) \
    .def("get_" #scatterer "_C" #name, &Experiment::get_##scatterer##_C##name, "Retrieves the scattering cross-section (C" #name ") for a" #scatterer)

#define COEFFICIENT_PROPERTY(scatterer, number) \
    .def("get_" #scatterer "_a" #number, &Experiment::get_##scatterer##_a##number, "Retrieves the electric multipole coefficient (a" #number ") coefficient for a " #scatterer) \
    .def("get_" #scatterer "_b" #number, &Experiment::get_##scatterer##_b##number, "Retrieves the magnetic multipole coefficient (b" #number ") coefficient for a " #scatterer)

#define COUPLING_AND_G(scatterer) \
     .def("get_" #scatterer "_g", &Experiment::get_##scatterer##_g, "Retrieves the asymmetry parameter (g) for a " #scatterer) \
     .def("get_" #scatterer "_coupling", &Experiment::get_##scatterer##_coupling, "Retrieves the coupling efficiency for a " #scatterer)


namespace py = pybind11;

PYBIND11_MODULE(Experiment, module) {
    module.doc() = "Interface for conducting Lorenz-Mie Theory (LMT) experiments within the PyMieSim package.";

    py::class_<Experiment>(module, "CppExperiment")
        .def(py::init<>(), "Constructs an Experiment object.")

        // Setup methods
        .def("set_detector", &Experiment::set_detector, "Configures the detector for the experiment.")
        .def("set_source", &Experiment::set_source, "Sets the light source for the experiment.")
        .def("set_sphere", &Experiment::set_sphere, "Defines a spherical scatterer for the experiment.")
        .def("set_cylinder", &Experiment::set_cylinder, "Defines a cylindrical scatterer for the experiment.")
        .def("set_coreshell", &Experiment::set_coreshell, "Defines a core-shell scatterer for the experiment.")


        // Sphere metrics
        EFFICIENCY_PROPERTY(sphere, sca)
        EFFICIENCY_PROPERTY(sphere, ext)
        EFFICIENCY_PROPERTY(sphere, abs)
        EFFICIENCY_PROPERTY(sphere, pr)
        EFFICIENCY_PROPERTY(sphere, forward)
        EFFICIENCY_PROPERTY(sphere, back)
        EFFICIENCY_PROPERTY(sphere, ratio)

        COEFFICIENT_PROPERTY(sphere, 1)
        COEFFICIENT_PROPERTY(sphere, 2)
        COEFFICIENT_PROPERTY(sphere, 3)

        COUPLING_AND_G(sphere)


        // Cylinder metrics
        EFFICIENCY_PROPERTY(cylinder, sca)
        EFFICIENCY_PROPERTY(cylinder, ext)
        EFFICIENCY_PROPERTY(cylinder, abs)

        COEFFICIENT_PROPERTY(cylinder, 11)
        COEFFICIENT_PROPERTY(cylinder, 21)
        COEFFICIENT_PROPERTY(cylinder, 12)
        COEFFICIENT_PROPERTY(cylinder, 22)
        COEFFICIENT_PROPERTY(cylinder, 13)
        COEFFICIENT_PROPERTY(cylinder, 23)

        COUPLING_AND_G(cylinder)


        // Coreshell metrics
        EFFICIENCY_PROPERTY(coreshell, sca)
        EFFICIENCY_PROPERTY(coreshell, ext)
        EFFICIENCY_PROPERTY(coreshell, abs)
        EFFICIENCY_PROPERTY(coreshell, pr)
        EFFICIENCY_PROPERTY(coreshell, forward)
        EFFICIENCY_PROPERTY(coreshell, back)
        EFFICIENCY_PROPERTY(coreshell, ratio)

        COEFFICIENT_PROPERTY(coreshell, 1)
        COEFFICIENT_PROPERTY(coreshell, 2)
        COEFFICIENT_PROPERTY(coreshell, 3)

        COUPLING_AND_G(coreshell)

        ;


}

