// #include <pybind11/pybind11.h>
// #include "sources.h"

// namespace py = pybind11;
// using namespace SOURCE;

// PYBIND11_MODULE(SourceInterface, module) {
//     module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

//     py::class_<BaseSource>(module, "BindedBaseSource")
//         .def(py::init<>())
//     ;

//     py::class_<Planewave, BaseSource>(module, "BindedPlanewave")
//         .def(
//             py::init<double, std::vector<complex128>, double>(),
//             py::arg("wavelength"),
//             py::arg("jones_vector"),
//             py::arg("amplitude"),
//             "Constructs a Planewave source with specified optical properties. ")

//         .def_readonly("wavelength", &Planewave::wavelength, "Wavelength of the source.")
//         .def_readonly("jones_vector", &Planewave::jones_vector, "Jones vector of the source.")
//         .def_readonly("amplitude", &Planewave::amplitude, "Electric field amplitude of the source.")
//         ;

//     py::class_<Gaussian, BaseSource>(module, "BindedGaussian")
//         .def(
//             py::init<double, std::vector<complex128>, double, double>(),
//             py::arg("wavelength"),
//             py::arg("jones_vector"),
//             py::arg("NA"),
//             py::arg("optical_power"),
//             "Constructs a Gaussian source with specified optical properties. ")

//         .def(py::init<>())

//         .def_readonly("wavelength", &Gaussian::wavelength, "Wavelength of the source.")
//         .def_readonly("jones_vector", &Gaussian::jones_vector, "Jones vector of the source.")
//         .def_readonly("amplitude", &Gaussian::amplitude, "Electric field amplitude of the source.")
//         ;
// }







#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <vector>
#include <complex>

using complex128 = std::complex<double>;

// Assuming SOURCE namespace contains your classes
namespace SOURCE {

    #define PI (double)3.14159265358979323846264338
    #define EPSILON0 (double)8.854187817620389e-12
    #define C (double)299792458.0

    class BaseSource {
    public:
        double wavelength;
        std::vector<complex128> jones_vector;
        double amplitude;
        double k;

        BaseSource() = default;
        BaseSource(double wavelength, std::vector<complex128> jones_vector, double amplitude)
            : wavelength(wavelength), jones_vector(jones_vector), amplitude(amplitude) {
            this->k = 2 * PI / this->wavelength;
        }
    };

    class Planewave : public BaseSource {
    public:
        Planewave() = default;
        Planewave(double wavelength, std::vector<complex128> jones_vector, double amplitude)
            : BaseSource(wavelength, jones_vector, amplitude) {}
    };

    class Gaussian : public BaseSource {
    public:
        double NA;
        double optical_power;

        Gaussian() = default;
        Gaussian(double wavelength, std::vector<complex128> jones_vector, double NA, double optical_power)
            : BaseSource(wavelength, jones_vector, 0.0), NA(NA), optical_power(optical_power) {
            }

    };

} // namespace SOURCE











namespace py = pybind11;
using namespace SOURCE;

PYBIND11_MODULE(SourceInterface, module) {
    module.doc() = "Lorenz-Mie Theory (LMT) C++ binding module for PyMieSim Python package.";

    py::class_<BaseSource>(module, "BindedBaseSource")
        .def(py::init<>());

    py::class_<Planewave, BaseSource>(module, "BindedPlanewave")
        .def(
            py::init<double, std::vector<complex128>, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("amplitude"),
            "Constructs a Planewave source with specified optical properties.")
        .def_readonly("wavelength", &Planewave::wavelength, "Wavelength of the source.")
        .def_readonly("jones_vector", &Planewave::jones_vector, "Jones vector of the source.")
        .def_readonly("amplitude", &Planewave::amplitude, "Electric field amplitude of the source.");

    py::class_<Gaussian, BaseSource>(module, "BindedGaussian")
        .def(
            py::init<double, std::vector<complex128>, double, double>(),
            py::arg("wavelength"),
            py::arg("jones_vector"),
            py::arg("NA"),
            py::arg("optical_power"),
            "Constructs a Gaussian source with specified optical properties.")
        .def(py::init<>())
        .def_readonly("wavelength", &Gaussian::wavelength, "Wavelength of the source.")
        .def_readonly("jones_vector", &Gaussian::jones_vector, "Jones vector of the source.")
        .def_readonly("amplitude", &Gaussian::amplitude, "Electric field amplitude of the source.")
        .def_readonly("NA", &Gaussian::NA, "Numerical Aperture of the source.")
        .def_readonly("optical_power", &Gaussian::optical_power, "Optical power of the source.");
}
