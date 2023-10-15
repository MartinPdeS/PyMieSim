#include <pybind11/pybind11.h>

#include "sets.cpp"


PYBIND11_MODULE(Sets, module)
{
      pybind11::class_<SPHERE::Set>(module, "CppSphereSet")
      .def(pybind11::init<DVector&, std::vector<complex128>&, DVector&>(),
           pybind11::arg("diameter"),
           pybind11::arg("index"),
           pybind11::arg("n_medium") )

      .def(pybind11::init<DVector&, std::vector<std::vector<complex128>>&, DVector&>(),
           pybind11::arg("diameter"),
           pybind11::arg("material"),
           pybind11::arg("n_medium") )
           ;

      pybind11::class_<CYLINDER::Set>(module, "CppCylinderSet")
      .def(pybind11::init<DVector&, CVector&, DVector&>(),
           pybind11::arg("diameter"),
           pybind11::arg("index"),
           pybind11::arg("n_medium")
           )

      .def(pybind11::init<DVector&, std::vector<std::vector<complex128>>&, DVector&>(),
           pybind11::arg("diameter"),
           pybind11::arg("material"),
           pybind11::arg("n_medium")
           );


    pybind11::class_<CORESHELL::Set>(module, "CppCoreShellSet")
    .def(pybind11::init<DVector&, DVector&, CVector&, CVector&, DVector&>(),
         pybind11::arg("core_diameter"),
         pybind11::arg("shell_width"),
         pybind11::arg("core_index"),
         pybind11::arg("shell_index"),
         pybind11::arg("n_medium") )

    .def(pybind11::init<DVector&, DVector&, CVector&, std::vector<CVector>&, DVector&>(),
         pybind11::arg("core_diameter"),
         pybind11::arg("shell_width"),
         pybind11::arg("core_index"),
         pybind11::arg("shell_material"),
         pybind11::arg("n_medium"))

    .def(pybind11::init<DVector&, DVector&, std::vector<CVector>&, CVector&, DVector&>(),
         pybind11::arg("core_diameter"),
         pybind11::arg("shell_width"),
         pybind11::arg("core_material"),
         pybind11::arg("shell_index"),
         pybind11::arg("n_medium") )

    .def(pybind11::init<DVector&, DVector&, std::vector<CVector>&, std::vector<CVector>&, DVector&>(),
         pybind11::arg("core_diameter"),
         pybind11::arg("shell_width"),
         pybind11::arg("core_material"),
         pybind11::arg("shell_material"),
         pybind11::arg("n_medium") )
         ;

      pybind11::class_<SOURCE::Set>(module, "CppSourceSet")
      .def(pybind11::init<DVector&, std::vector<CVector>&, DVector&>(),
           pybind11::arg("wavelength"),
           pybind11::arg("jones_vector"),
           pybind11::arg("amplitude") );


     pybind11::class_<DETECTOR::Set>(module, "CppDetectorSet")
     .def(pybind11::init<std::vector<CVector>&, DVector&, DVector&, DVector&, DVector&, bool&, bool&>(),
          pybind11::arg("scalarfield"),
          pybind11::arg("NA"),
          pybind11::arg("phi_offset"),
          pybind11::arg("gamma_offset"),
          pybind11::arg("polarization_filter"),
          pybind11::arg("coherent"),
          pybind11::arg("point_coupling")
          );

}






// -
