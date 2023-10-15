#include <pybind11/pybind11.h>
#include "sources.cpp"

PYBIND11_MODULE(DetectorInterface, module)
{
     module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";


      pybind11::class_<DETECTOR::Detector>(module, "BindedDetector")
      .def(pybind11::init<CVector&, double &, double &, double &, double &, bool &, bool &>(),
           pybind11::arg("ScalarField"),
           pybind11::arg("NA"),
           pybind11::arg("PhiOffset"),
           pybind11::arg("GammaOffset"),
           pybind11::arg("Filter"),
           pybind11::arg("Coherent"),
           pybind11::arg("PointCoupling")  //true = point ; false=mean
          )

      .def(pybind11::init<CVector &, double &, double &, double &, double &, bool &, bool &>(),
           pybind11::arg("ScalarField"),
           pybind11::arg("NA"),
           pybind11::arg("PhiOffset"),
           pybind11::arg("GammaOffset"),
           pybind11::arg("Filter"),
           pybind11::arg("Coherent"),
           pybind11::arg("PointCoupling")  //true = point ; false=mean
          )

       .def("CouplingSphere", &DETECTOR::Detector::Coupling<SPHERE::Scatterer>, pybind11::arg("Scatterer") )
       .def("CouplingCylinder", &DETECTOR::Detector::Coupling<CYLINDER::Scatterer>, pybind11::arg("Scatterer") )
       .def("CouplingCoreShell", &DETECTOR::Detector::Coupling<CORESHELL::Scatterer>, pybind11::arg("Scatterer") )
       ;

}







// -
