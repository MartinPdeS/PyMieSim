#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include "coordinates.h"
#include <utils/numpy_interface.h>
#include <pint/pint.h>


PYBIND11_MODULE(coordinates, module)
{
    py::object ureg = get_shared_ureg();


    pybind11::class_<Cartesian>(module, "Cartesian")
        .def(
            "to_spherical",
            &Cartesian::to_spherical,
            R"pbdoc(
                Converts the Cartesian coordinates to spherical coordinates.
                This method transforms the Cartesian coordinates of the mesh points into their corresponding spherical coordinates.
            )pbdoc"
        )
        .def_property_readonly(
            "x",
            [ureg](const Cartesian& self) -> pybind11::object {
                return pybind11::cast(self.x) * ureg("meter");
            },
            R"pbdoc(
                Returns x coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the x-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
        .def_property_readonly(
            "y",
            [ureg](const Cartesian& self) -> pybind11::object {
                return pybind11::cast(self.y) * ureg("meter");
            },
            R"pbdoc(
                Returns y coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the y-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
        .def_property_readonly(
            "z",
            [ureg](const Cartesian& self) -> pybind11::object {
                return pybind11::cast(self.z) * ureg("meter");
            },
            R"pbdoc(
                Returns z coordinates of points on the Cartesian mesh as a NumPy array.
                This property provides the z-coordinates of the mesh points in Cartesian coordinates.
            )pbdoc"
        )
    ;

    // ------------------ Bindings for SphericalCoordinate ------------------
    pybind11::class_<Spherical>(module, "Spherical")
        .def(
            "to_cartesian",
            &Spherical::to_cartesian,
            R"pbdoc(
                Converts the spherical coordinates to Cartesian coordinates.
                This method transforms the spherical coordinates of the mesh points into their corresponding Cartesian coordinates.
            )pbdoc"
        )
        .def_property_readonly(
            "r",
            [ureg](const Spherical& self) {
                const pybind11::ssize_t number_of_rows = self.shape[0];
                const pybind11::ssize_t number_of_columns = self.shape[1];

                return pybind11::array_t<double>(
                    {number_of_rows, number_of_columns},
                    {
                        static_cast<pybind11::ssize_t>(sizeof(double) * number_of_columns),
                        static_cast<pybind11::ssize_t>(sizeof(double))
                    },
                    self.r.data(),
                    pybind11::cast(self)
                ) * ureg("meter");
            },
            R"pbdoc(
                Returns radial distances of points on the spherical mesh as a NumPy array.
                This property provides the radial distances of the mesh points in spherical coordinates.
            )pbdoc"
        )
        .def_property_readonly(
            "phi",
            [ureg](const Spherical& self) {
                const pybind11::ssize_t number_of_rows = self.shape[0];
                const pybind11::ssize_t number_of_columns = self.shape[1];

                return pybind11::array_t<double>(
                    {number_of_rows, number_of_columns},
                    {
                        static_cast<pybind11::ssize_t>(sizeof(double) * number_of_columns),
                        static_cast<pybind11::ssize_t>(sizeof(double))
                    },
                    self.phi.data(),
                    pybind11::cast(self)
                ) * ureg("radian");
            },
            R"pbdoc(
                Returns azimuthal angles (phi) of points on the spherical mesh as a NumPy array.
                This property provides the azimuthal angles of the mesh points in spherical coordinates.
            )pbdoc"
        )
        .def_property_readonly(
            "theta",
            [ureg](const Spherical& self) {
                const pybind11::ssize_t number_of_rows = self.shape[0];
                const pybind11::ssize_t number_of_columns = self.shape[1];

                return pybind11::array_t<double>(
                    {number_of_rows, number_of_columns},
                    {
                        static_cast<pybind11::ssize_t>(sizeof(double) * number_of_columns),
                        static_cast<pybind11::ssize_t>(sizeof(double))
                    },
                    self.theta.data(),
                    pybind11::cast(self)
                ) * ureg("radian");
            },
            R"pbdoc(
                Returns polar angles (theta) of points on the spherical mesh as a NumPy array.
                This property provides the azimuthal angles of the mesh points in spherical coordinates.
            )pbdoc"
        )
    ;

    pybind11::class_<VectorField>(module, "VectorField")
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
