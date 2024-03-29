#include <pybind11/pybind11.h>

#include "experiment.cpp"



PYBIND11_MODULE(Experiment, module0)
{
  pybind11::class_<Experiment>(module0, "CppExperiment")
  .def(pybind11::init<>())

  .def("set_detector",               &Experiment::set_detector)
  .def("set_source",                 &Experiment::set_source)
  .def("set_sphere",                 &Experiment::set_sphere)
  .def("set_cylinder",               &Experiment::set_cylinder)
  .def("set_coreshell",              &Experiment::set_coreshell)

  .def("get_sphere_Qsca",            &Experiment::get_sphere_Qsca)
  .def("get_sphere_Qext",            &Experiment::get_sphere_Qext)
  .def("get_sphere_Qabs",            &Experiment::get_sphere_Qabs)
  .def("get_sphere_Qpr",             &Experiment::get_sphere_Qpr)
  .def("get_sphere_Qforward",        &Experiment::get_sphere_Qforward)
  .def("get_sphere_Qback",           &Experiment::get_sphere_Qback)

  .def("get_sphere_Csca",            &Experiment::get_sphere_Csca)
  .def("get_sphere_Cext",            &Experiment::get_sphere_Cext)
  .def("get_sphere_Cabs",            &Experiment::get_sphere_Cabs)
  .def("get_sphere_Cabs",            &Experiment::get_sphere_Cabs)
  .def("get_sphere_Cpr",             &Experiment::get_sphere_Cpr)
  .def("get_sphere_Cforward",        &Experiment::get_sphere_Cforward)
  .def("get_sphere_Cback",           &Experiment::get_sphere_Cback)

  .def("get_sphere_g",               &Experiment::get_sphere_g)
  .def("get_sphere_coupling",        &Experiment::get_sphere_coupling)

  .def("get_sphere_an",              &Experiment::get_sphere_an)
  .def("get_sphere_bn",              &Experiment::get_sphere_bn)
  .def("get_sphere_a1",              &Experiment::get_sphere_a1)
  .def("get_sphere_b1",              &Experiment::get_sphere_b1)
  .def("get_sphere_a2",              &Experiment::get_sphere_a2)
  .def("get_sphere_b2",              &Experiment::get_sphere_b2)
  .def("get_sphere_a3",              &Experiment::get_sphere_a3)
  .def("get_sphere_b3",              &Experiment::get_sphere_b3)

  .def("get_cylinder_Qsca",          &Experiment::get_cylinder_Qsca)
  .def("get_cylinder_Qext",          &Experiment::get_cylinder_Qext)
  .def("get_cylinder_Qabs",          &Experiment::get_cylinder_Qabs)
  .def("get_cylinder_Qpr",           &Experiment::get_cylinder_Qpr)
  .def("get_cylinder_Qforward",      &Experiment::get_cylinder_Qforward)
  .def("get_cylinder_Qback",         &Experiment::get_cylinder_Qback)

  .def("get_cylinder_Csca",          &Experiment::get_cylinder_Csca)
  .def("get_cylinder_Cext",          &Experiment::get_cylinder_Cext)
  .def("get_cylinder_Cabs",          &Experiment::get_cylinder_Cabs)
  .def("get_cylinder_Cabs",          &Experiment::get_cylinder_Cabs)
  .def("get_cylinder_Cpr",           &Experiment::get_cylinder_Cpr)
  .def("get_cylinder_Cforward",      &Experiment::get_cylinder_Cforward)
  .def("get_cylinder_Cback",         &Experiment::get_cylinder_Cback)

  .def("get_cylinder_g",             &Experiment::get_cylinder_g)
  .def("get_cylinder_coupling",      &Experiment::get_cylinder_coupling)

  .def("get_cylinder_a1n",           &Experiment::get_cylinder_a1n)
  .def("get_cylinder_a2n",           &Experiment::get_cylinder_a2n)
  .def("get_cylinder_b1n",           &Experiment::get_cylinder_b1n)
  .def("get_cylinder_b2n",           &Experiment::get_cylinder_b2n)
  .def("get_cylinder_a11",           &Experiment::get_cylinder_a11)
  .def("get_cylinder_a21",           &Experiment::get_cylinder_a21)
  .def("get_cylinder_b11",           &Experiment::get_cylinder_b11)
  .def("get_cylinder_b21",           &Experiment::get_cylinder_b21)
  .def("get_cylinder_a12",           &Experiment::get_cylinder_a12)
  .def("get_cylinder_a22",           &Experiment::get_cylinder_a22)
  .def("get_cylinder_b12",           &Experiment::get_cylinder_b12)
  .def("get_cylinder_b22",           &Experiment::get_cylinder_b22)
  .def("get_cylinder_a13",           &Experiment::get_cylinder_a13)
  .def("get_cylinder_a23",           &Experiment::get_cylinder_a23)
  .def("get_cylinder_b13",           &Experiment::get_cylinder_b13)
  .def("get_cylinder_b23",           &Experiment::get_cylinder_b23)

  .def("get_coreshell_Qsca",         &Experiment::get_coreshell_Qsca)
  .def("get_coreshell_Qext",         &Experiment::get_coreshell_Qext)
  .def("get_coreshell_Qabs",         &Experiment::get_coreshell_Qabs)
  .def("get_coreshell_Qpr",          &Experiment::get_coreshell_Qpr)
  .def("get_coreshell_Qforward",     &Experiment::get_coreshell_Qforward)
  .def("get_coreshell_Qback",        &Experiment::get_coreshell_Qback)

  .def("get_coreshell_Csca",         &Experiment::get_coreshell_Csca)
  .def("get_coreshell_Cext",         &Experiment::get_coreshell_Cext)
  .def("get_coreshell_Cabs",         &Experiment::get_coreshell_Cabs)
  .def("get_coreshell_Cabs",         &Experiment::get_coreshell_Cabs)
  .def("get_coreshell_Cpr",          &Experiment::get_coreshell_Cpr)
  .def("get_coreshell_Cforward",     &Experiment::get_coreshell_Cforward)
  .def("get_coreshell_Cback",        &Experiment::get_coreshell_Cback)

  .def("get_coreshell_g",            &Experiment::get_coreshell_g)
  .def("get_coreshell_coupling",     &Experiment::get_coreshell_coupling)

  .def("get_coreshell_an",           &Experiment::get_coreshell_an)
  .def("get_coreshell_bn",           &Experiment::get_coreshell_bn)
  .def("get_coreshell_a1",           &Experiment::get_coreshell_a1)
  .def("get_coreshell_b1",           &Experiment::get_coreshell_b1)
  .def("get_coreshell_a2",           &Experiment::get_coreshell_a2)
  .def("get_coreshell_b2",           &Experiment::get_coreshell_b2)
  .def("get_coreshell_a3",           &Experiment::get_coreshell_a3)
  .def("get_coreshell_b3",           &Experiment::get_coreshell_b3)

  ;

}
