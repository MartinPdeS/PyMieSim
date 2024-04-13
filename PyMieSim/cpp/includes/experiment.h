#pragma once

#include "definitions.cpp"
#include "numpy_interface.cpp"
#include "sources.cpp"
#include "sphere.cpp"
#include "cylinder.cpp"
#include "core_shell.cpp"
#include "detectors.cpp"


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

        size_t flatten_multi_index(const std::vector<size_t> &multi_index, const std::vector<size_t> &dimension) const;
        //--------------------------------------SPHERE------------------------------------
        template<typename Function>
        pybind11::array_t<complex128> get_sphere_coefficient(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_sphere_data(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_sphere_coefficient_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_sphere_coefficient_index(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_sphere_data_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_sphere_data_index(Function function, size_t max_order=0) const;

        pybind11::array_t<double> get_sphere_coupling() const;
        pybind11::array_t<double> get_sphere_coupling_material() const;
        pybind11::array_t<double> get_sphere_coupling_index() const;

        //--------------------------------------CYLINDER------------------------------------
        template<typename Function>
        pybind11::array_t<complex128> get_cylinder_coefficient(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_cylinder_data(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_cylinder_coefficient_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_cylinder_coefficient_index(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_cylinder_data_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_cylinder_data_index(Function function, size_t max_order=0) const;

        pybind11::array_t<double> get_cylinder_coupling() const;
        pybind11::array_t<double> get_cylinder_coupling_bound() const;
        pybind11::array_t<double> get_cylinder_coupling_unbound() const;

        //--------------------------------------CORESHELL------------------------------------
        template<typename Function>
        pybind11::array_t<complex128> get_coreshell_coefficient(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_coreshell_data(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_coreshell_coefficient_core_material_shell_index(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_coreshell_coefficient_core_index_shell_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_coreshell_coefficient_core_material_shell_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<complex128> get_coreshell_coefficient_core_index_shell_index(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_coreshell_data_core_material_shell_index(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_coreshell_data_core_index_shell_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_coreshell_data_core_material_shell_material(Function function, size_t max_order=0) const;

        template<typename Function>
        pybind11::array_t<double> get_coreshell_data_core_index_shell_index(Function function, size_t max_order=0) const;

        pybind11::array_t<double> get_coreshell_coupling() const;
        pybind11::array_t<double> get_coreshell_coupling_core_index_shell_index() const;
        pybind11::array_t<double> get_coreshell_coupling_core_material_shell_index() const;
        pybind11::array_t<double> get_coreshell_coupling_core_index_shell_material() const;
        pybind11::array_t<double> get_coreshell_coupling_core_material_shell_material() const;

        pybind11::array_t<double> get_sphere_Qsca() const { return get_sphere_data( &SPHERE::Scatterer::get_Qsca ) ; }
        pybind11::array_t<double> get_sphere_Qext() const { return get_sphere_data( &SPHERE::Scatterer::get_Qext ) ; }
        pybind11::array_t<double> get_sphere_Qabs() const { return get_sphere_data( &SPHERE::Scatterer::get_Qabs ) ; }
        pybind11::array_t<double> get_sphere_Qpr() const { return get_sphere_data( &SPHERE::Scatterer::get_Qpr ) ; }
        pybind11::array_t<double> get_sphere_Qback() const { return get_sphere_data( &SPHERE::Scatterer::get_Qback ) ; }
        pybind11::array_t<double> get_sphere_Qforward() const { return get_sphere_data( &SPHERE::Scatterer::get_Qforward ) ; }
        pybind11::array_t<double> get_sphere_Csca() const { return get_sphere_data( &SPHERE::Scatterer::get_Csca ) ; }
        pybind11::array_t<double> get_sphere_Cext() const { return get_sphere_data( &SPHERE::Scatterer::get_Cext ) ; }
        pybind11::array_t<double> get_sphere_Cabs() const { return get_sphere_data( &SPHERE::Scatterer::get_Cabs ) ; }
        pybind11::array_t<double> get_sphere_Cpr() const { return get_sphere_data( &SPHERE::Scatterer::get_Cpr ) ; }
        pybind11::array_t<double> get_sphere_Cback() const { return get_sphere_data( &SPHERE::Scatterer::get_Cback ) ; }
        pybind11::array_t<double> get_sphere_Cforward() const { return get_sphere_data( &SPHERE::Scatterer::get_Cforward ) ; }
        pybind11::array_t<double> get_sphere_g() const { return get_sphere_data( &SPHERE::Scatterer::get_g ) ; }

        pybind11::array_t<complex128> get_sphere_an(size_t max_order) const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, max_order ) ; }
        pybind11::array_t<complex128> get_sphere_bn(size_t max_order) const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, max_order ) ; }
        pybind11::array_t<complex128> get_sphere_a1() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 1 ) ; }
        pybind11::array_t<complex128> get_sphere_b1() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 1 ) ; }
        pybind11::array_t<complex128> get_sphere_a2() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 2 ) ; }
        pybind11::array_t<complex128> get_sphere_b2() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 2 ) ; }
        pybind11::array_t<complex128> get_sphere_a3() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_an, 3 ) ; }
        pybind11::array_t<complex128> get_sphere_b3() const { return get_sphere_coefficient( &SPHERE::Scatterer::get_bn, 3 ) ; }

        pybind11::array_t<double> get_cylinder_Qsca() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Qsca ) ; }
        pybind11::array_t<double> get_cylinder_Qext() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Qext ) ; }
        pybind11::array_t<double> get_cylinder_Qabs() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Qabs ) ; }
        pybind11::array_t<double> get_cylinder_Qpr() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Qpr ) ; }
        pybind11::array_t<double> get_cylinder_Qback() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Qback ) ; }
        pybind11::array_t<double> get_cylinder_Qforward() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Qforward ) ; }
        pybind11::array_t<double> get_cylinder_Csca() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Csca ) ; }
        pybind11::array_t<double> get_cylinder_Cext() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Cext ) ; }
        pybind11::array_t<double> get_cylinder_Cabs() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Cabs ) ; }
        pybind11::array_t<double> get_cylinder_Cpr() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Cpr ) ; }
        pybind11::array_t<double> get_cylinder_Cback() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Cback ) ; }
        pybind11::array_t<double> get_cylinder_Cforward() const { return get_cylinder_data( &CYLINDER::Scatterer::get_Cforward ) ; }
        pybind11::array_t<double> get_cylinder_g() const { return get_cylinder_data( &CYLINDER::Scatterer::get_g ) ; }

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

        pybind11::array_t<double> get_coreshell_Qsca() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Qsca ) ; }
        pybind11::array_t<double> get_coreshell_Qext() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Qext ) ; }
        pybind11::array_t<double> get_coreshell_Qabs() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Qabs ) ; }
        pybind11::array_t<double> get_coreshell_Qpr() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Qpr ) ; }
        pybind11::array_t<double> get_coreshell_Qback() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Qback ) ; }
        pybind11::array_t<double> get_coreshell_Qforward() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Qforward ) ; }
        pybind11::array_t<double> get_coreshell_Csca() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Csca ) ; }
        pybind11::array_t<double> get_coreshell_Cext() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Cext ) ; }
        pybind11::array_t<double> get_coreshell_Cabs() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Cabs ) ; }
        pybind11::array_t<double> get_coreshell_Cpr() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Cpr ) ; }
        pybind11::array_t<double> get_coreshell_Cback() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Cback ) ; }
        pybind11::array_t<double> get_coreshell_Cforward() const { return get_coreshell_data( &CORESHELL::Scatterer::get_Cforward ) ; }
        pybind11::array_t<double> get_coreshell_g() const { return get_coreshell_data( &CORESHELL::Scatterer::get_g ) ; }

        pybind11::array_t<complex128> get_coreshell_an(size_t max_order) const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, max_order ) ; }
        pybind11::array_t<complex128> get_coreshell_bn(size_t max_order) const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, max_order ) ; }
        pybind11::array_t<complex128> get_coreshell_a1() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 1 ) ; }
        pybind11::array_t<complex128> get_coreshell_b1() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 1 ) ; }
        pybind11::array_t<complex128> get_coreshell_a2() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 2 ) ; }
        pybind11::array_t<complex128> get_coreshell_b2() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 2 ) ; }
        pybind11::array_t<complex128> get_coreshell_a3() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_an, 3 ) ; }
        pybind11::array_t<complex128> get_coreshell_b3() const { return get_coreshell_coefficient( &CORESHELL::Scatterer::get_bn, 3 ) ; }
};


