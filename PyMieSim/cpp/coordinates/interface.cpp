#include <pybind11/pybind11.h>
#include "coordinates.h"
#include "../utils/numpy_interface.h"


void register_coordinates(pybind11::module_& module) {
    pybind11::class_<Cartesian>(module, "CARTESIANCOORDINATE")
        .def_property_readonly("x",
            [](const Cartesian& self) {
                return pybind11::array_t<double>(self.x.size(), self.x.data(), pybind11::cast(self));
            },
            R"pbdoc(
                Returns x coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the x-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
        .def_property_readonly("y",
            [](const Cartesian& self) {
                return pybind11::array_t<double>(self.y.size(), self.y.data(), pybind11::cast(self));
            },
            R"pbdoc(
                Returns y coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the y-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
        .def_property_readonly("z",
            [](const Cartesian& self) {
                return pybind11::array_t<double>(self.z.size(), self.z.data(), pybind11::cast(self));
            },
            R"pbdoc(
                Returns z coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the z-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
    ;

    // ------------------ Bindings for SphericalCoordinate ------------------
    pybind11::class_<Spherical>(module, "SPHERICALCOORDINATE")
        .def_property_readonly("r",
            [](const Spherical& self) {
                std::vector<size_t> shape = {self.r.size()};
                std::vector<size_t> strides = get_stride<double>(shape);

                return pybind11::array_t<double>(shape, strides, self.r.data(), pybind11::cast(self));
            },
            R"pbdoc(
                Returns radial distances of points on the spherical mesh as a NumPy array.
                This property provides the radial distances of the mesh points in spherical coordinates.
            )pbdoc"
        )
        .def_property_readonly("phi",
            [](const Spherical& self) {
                std::vector<size_t> shape = {self.phi.size()};
                std::vector<size_t> strides = get_stride<double>(shape);

                return pybind11::array_t<double>(shape, strides, self.phi.data(), pybind11::cast(self));
            },
            R"pbdoc(
                Returns azimuthal angles (phi) of points on the spherical mesh as a NumPy array.
                This property provides the azimuthal angles of the mesh points in spherical coordinates.
            )pbdoc"
        )
        .def_property_readonly("theta",
            [](const Spherical& self) {
                std::vector<size_t> shape = {self.theta.size()};
                std::vector<size_t> strides = get_stride<double>(shape);

                return pybind11::array_t<double>(shape, strides, self.theta.data(), pybind11::cast(self));
            },
            R"pbdoc(
                Returns polar angles (theta) of points on the spherical mesh as a NumPy array.
                This property provides the polar angles of the mesh points in spherical coordinates.
            )pbdoc"
        )
    ;

    pybind11::class_<VectorField>(module, "VECTORFIELD")
        .def_property_readonly("data",
             [](const VectorField& self) {
                std::vector<size_t> shape = {self.sampling, 3};
                return _vector_to_numpy(self.data, shape);
            },
            R"pbdoc(
                Returns the vector field data as a NumPy array.
                This property provides access to the raw vector field data.
            )pbdoc"
        )
        .def("_cpp_get_scalar_product",
            &VectorField::get_scalar_product,
            pybind11::arg("base"),
            R"pbdoc(
                Computes the scalar product with a base vector field.

                Parameters
                ----------
                base : VECTORFIELD
                    The base vector field to compute the scalar product with.

                Returns
                -------
                numpy.ndarray
                    A NumPy array containing the scalar products.
            )pbdoc"
        )
    ;
}
