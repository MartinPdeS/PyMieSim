#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "full_mesh.h"
#include <pint/pint.h>

namespace py = pybind11;

PYBIND11_MODULE(full_mesh, module)
{
    py::object ureg = get_shared_ureg();

    pybind11::class_<FullSteradian, std::shared_ptr<FullSteradian>>(module, "FullMesh")
    .def_readonly(
        "spherical",
        &FullSteradian::spherical,
        R"pbdoc(
            The spherical coordinates of the full mesh.

            This attribute contains the spherical coordinate representation of the full mesh, including the radial distance, polar angle, and azimuthal angle for each point in the mesh.
        )pbdoc"
    )
    .def_readonly(
        "cartesian",
        &FullSteradian::cartesian,
        R"pbdoc(
            The Cartesian coordinates of the full mesh.
            This attribute contains the Cartesian coordinate representation of the full mesh, including the x, y, and z coordinates for each point in the mesh.

        )pbdoc"
    )
    .def_readonly(
        "spherical_mesh",
        &FullSteradian::spherical_mesh,
        R"pbdoc(
            The spherical coordinates of the full mesh.
            This attribute contains the spherical coordinate representation of the full mesh, including the radial distance, polar angle, and azimuthal angle for each point in the mesh.
        )pbdoc"
    )
    .def_readonly(
        "cartesian_mesh",
        &FullSteradian::cartesian_mesh,
        R"pbdoc(
            The Cartesian coordinates of the full mesh.
            This attribute contains the Cartesian coordinate representation of the full mesh, including the x, y, and z coordinates for each point in the mesh.
        )pbdoc"
    );


}
