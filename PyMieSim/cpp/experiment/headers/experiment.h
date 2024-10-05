#pragma once


#include "utils/numpy_interface.cpp"
#include "single/includes/sources.cpp"
#include "single/headers/core_shell.h"
#include "single/headers/sphere.h"
#include "single/headers/cylinder.h"
#include "single/includes/detectors.cpp"


typedef std::complex<double> complex128;

#define DEFINE_SPHERE_FUNCTION(dtype, name) \
    pybind11::array_t<dtype> get_sphere_##name() const { return get_sphere_data<dtype>(&SPHERE::Scatterer::get_##name); }

#define DEFINE_CYLINDER_FUNCTION(dtype, name) \
    pybind11::array_t<dtype> get_cylinder_##name() const { return get_cylinder_data<dtype>(&CYLINDER::Scatterer::get_##name); }

#define DEFINE_CORESHELL_FUNCTION(dtype, name) \
    pybind11::array_t<dtype> get_coreshell_##name() const { return get_coreshell_data<dtype>(&CORESHELL::Scatterer::get_##name); }


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

    static size_t flatten_multi_index(const std::vector<size_t>& multi_index, const std::vector<size_t>& dimensions) {
        size_t flatten_index = 0;
        size_t stride = 1;

        const size_t* multi_index_ptr = multi_index.data();
        const size_t* dimensions_ptr = dimensions.data();

        // Iterate from the last dimension to the first
        for (int i = dimensions.size() - 1; i >= 0; --i) {
            flatten_index += multi_index_ptr[i] * stride;
            stride *= dimensions_ptr[i];
        }

        return flatten_index;
    }

        //--------------------------------------SPHERE------------------------------------
        template<typename dtype, typename Function> pybind11::array_t<dtype> get_sphere_data(Function function) const;

        pybind11::array_t<double> get_sphere_coupling() const;

        DEFINE_SPHERE_FUNCTION(complex128, a1)
        DEFINE_SPHERE_FUNCTION(complex128, a2)
        DEFINE_SPHERE_FUNCTION(complex128, a3)
        DEFINE_SPHERE_FUNCTION(complex128, b1)
        DEFINE_SPHERE_FUNCTION(complex128, b2)
        DEFINE_SPHERE_FUNCTION(complex128, b3)

        DEFINE_SPHERE_FUNCTION(double, Qsca)
        DEFINE_SPHERE_FUNCTION(double, Qext)
        DEFINE_SPHERE_FUNCTION(double, Qabs)
        DEFINE_SPHERE_FUNCTION(double, Qpr)
        DEFINE_SPHERE_FUNCTION(double, Qback)
        DEFINE_SPHERE_FUNCTION(double, Qforward)
        DEFINE_SPHERE_FUNCTION(double, Qratio)
        DEFINE_SPHERE_FUNCTION(double, Csca)
        DEFINE_SPHERE_FUNCTION(double, Cext)
        DEFINE_SPHERE_FUNCTION(double, Cabs)
        DEFINE_SPHERE_FUNCTION(double, Cpr)
        DEFINE_SPHERE_FUNCTION(double, Cback)
        DEFINE_SPHERE_FUNCTION(double, Cratio)
        DEFINE_SPHERE_FUNCTION(double, Cforward)
        DEFINE_SPHERE_FUNCTION(double, g)

        //--------------------------------------CYLINDER------------------------------------
        template<typename dtype, typename Function> pybind11::array_t<dtype> get_cylinder_data(Function function) const;

        pybind11::array_t<double> get_cylinder_coupling() const;

        DEFINE_CYLINDER_FUNCTION(complex128, a11)
        DEFINE_CYLINDER_FUNCTION(complex128, a12)
        DEFINE_CYLINDER_FUNCTION(complex128, a13)
        DEFINE_CYLINDER_FUNCTION(complex128, a21)
        DEFINE_CYLINDER_FUNCTION(complex128, a22)
        DEFINE_CYLINDER_FUNCTION(complex128, a23)
        DEFINE_CYLINDER_FUNCTION(complex128, b11)
        DEFINE_CYLINDER_FUNCTION(complex128, b12)
        DEFINE_CYLINDER_FUNCTION(complex128, b13)
        DEFINE_CYLINDER_FUNCTION(complex128, b21)
        DEFINE_CYLINDER_FUNCTION(complex128, b22)
        DEFINE_CYLINDER_FUNCTION(complex128, b23)

        DEFINE_CYLINDER_FUNCTION(double, Qsca)
        DEFINE_CYLINDER_FUNCTION(double, Qext)
        DEFINE_CYLINDER_FUNCTION(double, Qabs)
        DEFINE_CYLINDER_FUNCTION(double, Csca)
        DEFINE_CYLINDER_FUNCTION(double, Cext)
        DEFINE_CYLINDER_FUNCTION(double, Cabs)
        DEFINE_CYLINDER_FUNCTION(double, g)

        //--------------------------------------CORESHELL------------------------------------
        template<typename dtype, typename Function> pybind11::array_t<dtype> get_coreshell_data(Function function) const;

        pybind11::array_t<double> get_coreshell_coupling() const;

        DEFINE_CORESHELL_FUNCTION(complex128, a1)
        DEFINE_CORESHELL_FUNCTION(complex128, a2)
        DEFINE_CORESHELL_FUNCTION(complex128, a3)
        DEFINE_CORESHELL_FUNCTION(complex128, b1)
        DEFINE_CORESHELL_FUNCTION(complex128, b2)
        DEFINE_CORESHELL_FUNCTION(complex128, b3)

        DEFINE_CORESHELL_FUNCTION(double, Qsca)
        DEFINE_CORESHELL_FUNCTION(double, Qext)
        DEFINE_CORESHELL_FUNCTION(double, Qabs)
        DEFINE_CORESHELL_FUNCTION(double, Qpr)
        DEFINE_CORESHELL_FUNCTION(double, Qback)
        DEFINE_CORESHELL_FUNCTION(double, Qforward)
        DEFINE_CORESHELL_FUNCTION(double, Qratio)
        DEFINE_CORESHELL_FUNCTION(double, Csca)
        DEFINE_CORESHELL_FUNCTION(double, Cext)
        DEFINE_CORESHELL_FUNCTION(double, Cabs)
        DEFINE_CORESHELL_FUNCTION(double, Cpr)
        DEFINE_CORESHELL_FUNCTION(double, Cback)
        DEFINE_CORESHELL_FUNCTION(double, Cratio)
        DEFINE_CORESHELL_FUNCTION(double, Cforward)
        DEFINE_CORESHELL_FUNCTION(double, g)

};


