#include <pybind11/pybind11.h>
#include "experiment.cpp" // Make sure this includes the definition of the Experiment class

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
        // Downward are the sphere efficiencies
        .def("get_sphere_Qsca", &Experiment::get_sphere_Qsca, "Retrieves the scattering efficiency (Qsca) for a sphere.")
        .def("get_sphere_Qext", &Experiment::get_sphere_Qext, "Retrieves the extinction efficiency (Qext) for a sphere.")
        .def("get_sphere_Qabs", &Experiment::get_sphere_Qabs, "Retrieves the absorption efficiency (Qabs) for a sphere.")
        .def("get_sphere_Qpr", &Experiment::get_sphere_Qpr, "Retrieves the radiation pressure efficiency (Qpr) for a sphere.")
        .def("get_sphere_Qforward", &Experiment::get_sphere_Qforward, "Retrieves the forward scattering efficiency (Qforward) for a sphere.")
        .def("get_sphere_Qback", &Experiment::get_sphere_Qback, "Retrieves the backscattering efficiency (Qback) for a sphere.")
        .def("get_sphere_Qratio", &Experiment::get_sphere_Qratio, "Retrieves the ratio between forward and backward efficiencies for a sphere.")
        // Downward are the sphere cross-sections
        .def("get_sphere_Csca", &Experiment::get_sphere_Csca, "Retrieves the scattering cross-section (Csca) for a sphere.")
        .def("get_sphere_Cext", &Experiment::get_sphere_Cext, "Retrieves the extinction cross-section (Cext) for a sphere.")
        .def("get_sphere_Cabs", &Experiment::get_sphere_Cabs, "Retrieves the absorption cross-section (Cabs) for a sphere.")
        .def("get_sphere_Cpr", &Experiment::get_sphere_Cpr, "Retrieves the radiation pressure cross-section (Cpr) for a sphere.")
        .def("get_sphere_Cforward", &Experiment::get_sphere_Cforward, "Retrieves the forward scattering cross-section (Cforward) for a sphere.")
        .def("get_sphere_Cback", &Experiment::get_sphere_Cback, "Retrieves the backscattering cross-section (Cback) for a sphere.")
        .def("get_sphere_Cratio", &Experiment::get_sphere_Cratio, "Retrieves the ratio between forward and backward cross-section for a sphere.")
        // Downward are the sphere extra parameters
        .def("get_sphere_g", &Experiment::get_sphere_g, "Retrieves the asymmetry parameter (g) for a sphere.")
        .def("get_sphere_coupling", &Experiment::get_sphere_coupling, "Retrieves the coupling efficiency for a sphere.")
        // Sphere coefficient retrievals
        .def("get_sphere_an", &Experiment::get_sphere_an, "Retrieves the an coefficient for a sphere.")
        .def("get_sphere_bn", &Experiment::get_sphere_bn, "Retrieves the bn coefficient for a sphere.")
        .def("get_sphere_a1", &Experiment::get_sphere_a1, "Retrieves the a1 coefficient for a sphere.")
        .def("get_sphere_b1", &Experiment::get_sphere_b1, "Retrieves the b1 coefficient for a sphere.")
        .def("get_sphere_a2", &Experiment::get_sphere_a2, "Retrieves the a2 coefficient for a sphere.")
        .def("get_sphere_b2", &Experiment::get_sphere_b2, "Retrieves the b2 coefficient for a sphere.")
        .def("get_sphere_a3", &Experiment::get_sphere_a3, "Retrieves the a3 coefficient for a sphere.")
        .def("get_sphere_b3", &Experiment::get_sphere_b3, "Retrieves the b3 coefficient for a sphere.")

        // Cylinder metrics
        // Downward are the cylinder efficiencies
        .def("get_cylinder_Qsca", &Experiment::get_cylinder_Qsca, "Retrieves the scattering efficiency (Qsca) for a cylinder.")
        .def("get_cylinder_Qext", &Experiment::get_cylinder_Qext, "Retrieves the extinction efficiency (Qext) for a cylinder.")
        .def("get_cylinder_Qabs", &Experiment::get_cylinder_Qabs, "Retrieves the absorption efficiency (Qabs) for a cylinder.")
        // Downward are the cylinder cross-sections
        .def("get_cylinder_Csca", &Experiment::get_cylinder_Csca, "Retrieves the scattering cross-section (Csca) for a cylinder.")
        .def("get_cylinder_Cext", &Experiment::get_cylinder_Cext, "Retrieves the extinction cross-section (Cext) for a cylinder.")
        .def("get_cylinder_Cabs", &Experiment::get_cylinder_Cabs, "Retrieves the absorption cross-section (Cabs) for a cylinder.")
        // Downward are the cylinder extra parameters
        .def("get_cylinder_g", &Experiment::get_cylinder_g, "Retrieves the asymmetry parameter (g) for a cylinder.")
        .def("get_cylinder_coupling", &Experiment::get_cylinder_coupling, "Retrieves the coupling efficiency for a cylinder.")

        // Cylinder coefficient retrievals
        .def("get_cylinder_a1n", &Experiment::get_cylinder_a1n, "Retrieves the a1n coefficient for a cylinder.")
        .def("get_cylinder_a2n", &Experiment::get_cylinder_a2n, "Retrieves the a2n coefficient for a cylinder.")
        .def("get_cylinder_b1n", &Experiment::get_cylinder_b1n, "Retrieves the b1n coefficient for a cylinder.")
        .def("get_cylinder_b2n", &Experiment::get_cylinder_b2n, "Retrieves the b2n coefficient for a cylinder.")
        .def("get_cylinder_a11", &Experiment::get_cylinder_a11, "Retrieves the a11 coefficient for a cylinder.")
        .def("get_cylinder_a21", &Experiment::get_cylinder_a21, "Retrieves the a21 coefficient for a cylinder.")
        .def("get_cylinder_b11", &Experiment::get_cylinder_b11, "Retrieves the b11 coefficient for a cylinder.")
        .def("get_cylinder_b21", &Experiment::get_cylinder_b21, "Retrieves the b21 coefficient for a cylinder.")
        .def("get_cylinder_a12", &Experiment::get_cylinder_a12, "Retrieves the a12 coefficient for a cylinder.")
        .def("get_cylinder_a22", &Experiment::get_cylinder_a22, "Retrieves the a22 coefficient for a cylinder.")
        .def("get_cylinder_b12", &Experiment::get_cylinder_b12, "Retrieves the b12 coefficient for a cylinder.")
        .def("get_cylinder_b22", &Experiment::get_cylinder_b22, "Retrieves the b22 coefficient for a cylinder.")
        .def("get_cylinder_a13", &Experiment::get_cylinder_a13, "Retrieves the a13 coefficient for a cylinder.")
        .def("get_cylinder_a23", &Experiment::get_cylinder_a23, "Retrieves the a23 coefficient for a cylinder.")
        .def("get_cylinder_b13", &Experiment::get_cylinder_b13, "Retrieves the b13 coefficient for a cylinder.")
        .def("get_cylinder_b23", &Experiment::get_cylinder_b23, "Retrieves the b23 coefficient for a cylinder.")

        // Coreshell metrics
        // Downward are the core/shell efficiencies
        .def("get_coreshell_Qsca", &Experiment::get_coreshell_Qsca, "Retrieves the scattering efficiency (Qsca) for a coreshell.")
        .def("get_coreshell_Qext", &Experiment::get_coreshell_Qext, "Retrieves the extinction efficiency (Qext) for a coreshell.")
        .def("get_coreshell_Qabs", &Experiment::get_coreshell_Qabs, "Retrieves the absorption efficiency (Qabs) for a coreshell.")
        .def("get_coreshell_Qpr", &Experiment::get_coreshell_Qpr, "Retrieves the radiation pressure efficiency (Qpr) for a coreshell.")
        .def("get_coreshell_Qforward", &Experiment::get_coreshell_Qforward, "Retrieves the forward scattering efficiency (Qforward) for a coreshell.")
        .def("get_coreshell_Qback", &Experiment::get_coreshell_Qback, "Retrieves the backscattering efficiency (Qback) for a coreshell.")
        .def("get_coreshell_Qratio", &Experiment::get_coreshell_Qratio, "Retrieves the ratio between forward and backward efficiencies for a coreshell.")
        // Downward are the core/shell cross-sections
        .def("get_coreshell_Csca", &Experiment::get_coreshell_Csca, "Retrieves the scattering cross-section (Csca) for a coreshell.")
        .def("get_coreshell_Cext", &Experiment::get_coreshell_Cext, "Retrieves the extinction cross-section (Cext) for a coreshell.")
        .def("get_coreshell_Cabs", &Experiment::get_coreshell_Cabs, "Retrieves the absorption cross-section (Cabs) for a coreshell.")
        .def("get_coreshell_Cpr", &Experiment::get_coreshell_Cpr, "Retrieves the radiation pressure cross-section (Cpr) for a coreshell.")
        .def("get_coreshell_Cforward", &Experiment::get_coreshell_Cforward, "Retrieves the forward scattering cross-section (Cforward) for a coreshell.")
        .def("get_coreshell_Cback", &Experiment::get_coreshell_Cback, "Retrieves the backscattering cross-section (Cback) for a coreshell.")
        .def("get_coreshell_Cratio", &Experiment::get_coreshell_Cratio, "Retrieves the ratio between forward and backward cross-section for a coreshell.")
        // Downward are the core/shell extra parameters
        .def("get_coreshell_g", &Experiment::get_coreshell_g, "Retrieves the asymmetry parameter (g) for a coreshell.")
        .def("get_coreshell_coupling", &Experiment::get_coreshell_coupling, "Retrieves the coupling efficiency for a coreshell.")

        // Coreshell coefficient retrievals
        .def("get_coreshell_an", &Experiment::get_coreshell_an, "Retrieves the an coefficient for a coreshell.")
        .def("get_coreshell_bn", &Experiment::get_coreshell_bn, "Retrieves the bn coefficient for a coreshell.")
        .def("get_coreshell_a1", &Experiment::get_coreshell_a1, "Retrieves the a1 coefficient for a coreshell.")
        .def("get_coreshell_b1", &Experiment::get_coreshell_b1, "Retrieves the b1 coefficient for a coreshell.")
        .def("get_coreshell_a2", &Experiment::get_coreshell_a2, "Retrieves the a2 coefficient for a coreshell.")
        .def("get_coreshell_b2", &Experiment::get_coreshell_b2, "Retrieves the b2 coefficient for a coreshell.")
        .def("get_coreshell_a3", &Experiment::get_coreshell_a3, "Retrieves the a3 coefficient for a coreshell.")
        .def("get_coreshell_b3", &Experiment::get_coreshell_b3, "Retrieves the b3 coefficient for a coreshell.");

}

