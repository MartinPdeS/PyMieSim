#include <pybind11/pybind11.h>
#include "detectors.cpp"
#include "cylinder.cpp"
#include "sphere.cpp"
#include "core_shell.cpp"

PYBIND11_MODULE(DetectorInterface, module)
{
     module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";


      pybind11::class_<DETECTOR::Detector>(module, "BindedDetector")
      .def(pybind11::init<CVector&, double &, double &, double &, double &, bool &, bool &>(),
           pybind11::arg("scalar_field"),
           pybind11::arg("NA"),
           pybind11::arg("phi_offset"),
           pybind11::arg("gamma_offset"),
           pybind11::arg("polarization_filter"),
           pybind11::arg("coherent"),
           pybind11::arg("point_coupling")  //true = point ; false=mean
          )

      .def(pybind11::init<CVector &, double &, double &, double &, double &, bool &, bool &>(),
           pybind11::arg("scalar_field"),
           pybind11::arg("NA"),
           pybind11::arg("phi_offset"),
           pybind11::arg("gamma_offset"),
           pybind11::arg("polarization_filter"),
           pybind11::arg("coherent"),
           pybind11::arg("point_coupling")  //true = point ; false=mean
          )

       .def("CouplingSphere", &DETECTOR::Detector::Coupling<SPHERE::Scatterer>, pybind11::arg("Scatterer") )
       .def("CouplingCylinder", &DETECTOR::Detector::Coupling<CYLINDER::Scatterer>, pybind11::arg("Scatterer") )
       .def("CouplingCoreShell", &DETECTOR::Detector::Coupling<CORESHELL::Scatterer>, pybind11::arg("Scatterer") )
       ;

}







// -
