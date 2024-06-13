#pragma once

#include "definitions.cpp"
#include "numpy_interface.cpp"
#include "sources.cpp"
#include "sphere.cpp"
#include "cylinder.cpp"
#include "core_shell.cpp"
#include "detectors.cpp"

#define DEFINE_SPHERE_FUNCTION(name) \
    pybind11::array_t<double> get_sphere_##name() const { return get_sphere_data(&SPHERE::Scatterer::get_##name); }

#define DEFINE_CYLINDER_FUNCTION(name) \
    pybind11::array_t<double> get_cylinder_##name() const { return get_cylinder_data(&CYLINDER::Scatterer::get_##name); }

#define DEFINE_CORESHELL_FUNCTION(name) \
    pybind11::array_t<double> get_coreshell_##name() const { return get_coreshell_data(&CORESHELL::Scatterer::get_##name); }


class Experiment
{
    public:
        SPHERE::Set sphereSet;
        CYLINDER::Set cylinderSet;
        CORESHELL::Set coreshellSet;
        DETECTOR::Set detectorSet;
        SOURCE::Set sourceSet;

        Experiment() = default;

        void set_sphere(SPHERE::Set& set) { sphereSet = set; }
        void set_cylinder(CYLINDER::Set& set) { cylinderSet = set; }
        void set_coreshell(CORESHELL::Set& set) { coreshellSet = set; }
        void set_source(SOURCE::Set &set) { sourceSet = set; }
        void set_detector(DETECTOR::Set &set) { detectorSet = set; }

        static size_t flatten_multi_index(const std::vector<size_t>& multi_index, const std::vector<size_t>& dimensions) { // Trust chatGPT on that one
            size_t flatten_index = 0;
            size_t stride = 1;

            // Iterate from the last dimension to the first
            for (int i = dimensions.size() - 1; i >= 0; --i) {
                flatten_index += multi_index[i] * stride;
                stride *= dimensions[i];
            }

            return flatten_index;
        }

        //--------------------------------------SPHERE------------------------------------
        template<typename Function> pybind11::array_t<complex128> get_sphere_coefficient(Function function, size_t max_order=0) const;

        template<typename Function> pybind11::array_t<double> get_sphere_data(Function function) const;

        pybind11::array_t<double> get_sphere_coupling() const;

        DEFINE_SPHERE_FUNCTION(Qsca)
        DEFINE_SPHERE_FUNCTION(Qext)
        DEFINE_SPHERE_FUNCTION(Qabs)
        DEFINE_SPHERE_FUNCTION(Qpr)
        DEFINE_SPHERE_FUNCTION(Qback)
        DEFINE_SPHERE_FUNCTION(Qforward)
        DEFINE_SPHERE_FUNCTION(Qratio)
        DEFINE_SPHERE_FUNCTION(Csca)
        DEFINE_SPHERE_FUNCTION(Cext)
        DEFINE_SPHERE_FUNCTION(Cabs)
        DEFINE_SPHERE_FUNCTION(Cpr)
        DEFINE_SPHERE_FUNCTION(Cback)
        DEFINE_SPHERE_FUNCTION(Cratio)
        DEFINE_SPHERE_FUNCTION(Cforward)
        DEFINE_SPHERE_FUNCTION(g)

        //--------------------------------------CYLINDER------------------------------------
        template<typename Function> pybind11::array_t<complex128> get_cylinder_coefficient(Function function, size_t max_order=0) const;

        template<typename Function> pybind11::array_t<double> get_cylinder_data(Function function) const;

        pybind11::array_t<double> get_cylinder_coupling() const;

        DEFINE_CYLINDER_FUNCTION(Qsca)
        DEFINE_CYLINDER_FUNCTION(Qext)
        DEFINE_CYLINDER_FUNCTION(Qabs)
        DEFINE_CYLINDER_FUNCTION(Qpr)
        DEFINE_CYLINDER_FUNCTION(Qback)
        DEFINE_CYLINDER_FUNCTION(Qforward)
        DEFINE_CYLINDER_FUNCTION(Csca)
        DEFINE_CYLINDER_FUNCTION(Cext)
        DEFINE_CYLINDER_FUNCTION(Cabs)
        DEFINE_CYLINDER_FUNCTION(Cpr)
        DEFINE_CYLINDER_FUNCTION(Cback)
        DEFINE_CYLINDER_FUNCTION(Cforward)
        DEFINE_CYLINDER_FUNCTION(g)

        //--------------------------------------CORESHELL------------------------------------
        template<typename Function> pybind11::array_t<complex128> get_coreshell_coefficient(Function function, size_t max_order=0) const;

        template<typename Function> pybind11::array_t<double> get_coreshell_data(Function function) const;

        pybind11::array_t<double> get_coreshell_coupling() const;

        DEFINE_CORESHELL_FUNCTION(Qsca)
        DEFINE_CORESHELL_FUNCTION(Qext)
        DEFINE_CORESHELL_FUNCTION(Qabs)
        DEFINE_CORESHELL_FUNCTION(Qpr)
        DEFINE_CORESHELL_FUNCTION(Qback)
        DEFINE_CORESHELL_FUNCTION(Qforward)
        DEFINE_CORESHELL_FUNCTION(Qratio)
        DEFINE_CORESHELL_FUNCTION(Csca)
        DEFINE_CORESHELL_FUNCTION(Cext)
        DEFINE_CORESHELL_FUNCTION(Cabs)
        DEFINE_CORESHELL_FUNCTION(Cpr)
        DEFINE_CORESHELL_FUNCTION(Cback)
        DEFINE_CORESHELL_FUNCTION(Cratio)
        DEFINE_CORESHELL_FUNCTION(Cforward)
        DEFINE_CORESHELL_FUNCTION(g)


        pybind11::array_t<complex128> get_sphere_an(size_t max_order) const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, max_order ) ; }
        pybind11::array_t<complex128> get_sphere_bn(size_t max_order) const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, max_order ) ; }
        pybind11::array_t<complex128> get_sphere_a1() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 1 ) ; }
        pybind11::array_t<complex128> get_sphere_b1() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 1 ) ; }
        pybind11::array_t<complex128> get_sphere_a2() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 2 ) ; }
        pybind11::array_t<complex128> get_sphere_b2() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 2 ) ; }
        pybind11::array_t<complex128> get_sphere_a3() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 3 ) ; }
        pybind11::array_t<complex128> get_sphere_b3() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 3 ) ; }

        pybind11::array_t<complex128> get_cylinder_a1n(size_t max_order) const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, max_order ) ; }
        pybind11::array_t<complex128> get_cylinder_b1n(size_t max_order) const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, max_order ) ; }
        pybind11::array_t<complex128> get_cylinder_a2n(size_t max_order) const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, max_order ) ; }
        pybind11::array_t<complex128> get_cylinder_b2n(size_t max_order) const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, max_order ) ; }
        pybind11::array_t<complex128> get_cylinder_a11() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, 1 ) ; }
        pybind11::array_t<complex128> get_cylinder_b11() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, 1 ) ; }
        pybind11::array_t<complex128> get_cylinder_a21() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, 1 ) ; }
        pybind11::array_t<complex128> get_cylinder_b21() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, 1 ) ; }
        pybind11::array_t<complex128> get_cylinder_a12() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, 2 ) ; }
        pybind11::array_t<complex128> get_cylinder_b12() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, 2 ) ; }
        pybind11::array_t<complex128> get_cylinder_a22() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, 2 ) ; }
        pybind11::array_t<complex128> get_cylinder_b22() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, 2 ) ; }
        pybind11::array_t<complex128> get_cylinder_a13() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a1n, 3 ) ; }
        pybind11::array_t<complex128> get_cylinder_b13() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b1n, 3 ) ; }
        pybind11::array_t<complex128> get_cylinder_a23() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_a2n, 3 ) ; }
        pybind11::array_t<complex128> get_cylinder_b23() const { return get_cylinder_coefficient( &CYLINDER::Scatterer::get_b2n, 3 ) ; }

        pybind11::array_t<complex128> get_coreshell_an(size_t max_order) const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, max_order ) ; }
        pybind11::array_t<complex128> get_coreshell_bn(size_t max_order) const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, max_order ) ; }
        pybind11::array_t<complex128> get_coreshell_a1() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 1 ) ; }
        pybind11::array_t<complex128> get_coreshell_b1() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 1 ) ; }
        pybind11::array_t<complex128> get_coreshell_a2() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 2 ) ; }
        pybind11::array_t<complex128> get_coreshell_b2() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 2 ) ; }
        pybind11::array_t<complex128> get_coreshell_a3() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 3 ) ; }
        pybind11::array_t<complex128> get_coreshell_b3() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 3 ) ; }
};


