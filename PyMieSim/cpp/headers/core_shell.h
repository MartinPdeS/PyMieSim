#pragma once

#include "base_scatterer.cpp"
#include "numpy_interface.cpp"
#include "VSH.cpp"

  namespace CORESHELL
  {

    class Set
    {
        public:
            std::vector<double> core_diameter;
            std::vector<double> shell_width;
            std::vector<double> n_medium;
            std::vector<complex128> core_index;
            std::vector<complex128> shell_index;
            std::vector<std::vector<complex128>> core_material;
            std::vector<std::vector<complex128>> shell_material;
            bool core_is_material;
            bool shell_is_material;

            std::vector<size_t> shape;

            Set() = default;

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<complex128> &core_index,
                const std::vector<complex128> &shell_index,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter), shell_width(shell_width), core_index(core_index),
                shell_index(shell_index), n_medium(n_medium), core_is_material(false), shell_is_material(false)
            {
                this->shape = {
                    this->core_diameter.size(),
                    this->shell_width.size(),
                    this->core_index.size(),
                    this->shell_index.size(),
                    this->n_medium.size()
                };
            }

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<complex128> &core_index,
                const std::vector<std::vector<complex128>> &shell_material,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter), shell_width(shell_width), core_index(core_index),
                shell_material(shell_material), n_medium(n_medium), core_is_material(false), shell_is_material(true)
            {
                this->shape = {
                    this->core_diameter.size(),
                    this->shell_width.size(),
                    this->core_index.size(),
                    this->shell_material.size(),
                    this->n_medium.size()
                };
            }

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<std::vector<complex128>> &core_material,
                const std::vector<complex128> &shell_index,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter), shell_width(shell_width), shell_index(shell_index),
                core_material(core_material), n_medium(n_medium), core_is_material(true), shell_is_material(false)
            {
                this->shape = {
                    this->core_diameter.size(),
                    this->shell_width.size(),
                    this->core_material.size(),
                    this->shell_index.size(),
                    this->n_medium.size()
                };
            }

            Set(
                const std::vector<double> &core_diameter,
                const std::vector<double> &shell_width,
                const std::vector<std::vector<complex128>> &core_material,
                const std::vector<std::vector<complex128>> &shell_material,
                const std::vector<double> &n_medium)
            :
                core_diameter(core_diameter), shell_width(shell_width), core_material(core_material),
                shell_material(shell_material), n_medium(n_medium), core_is_material(true), shell_is_material(true)
            {
                this->shape = {
                    this->core_diameter.size(),
                    this->shell_width.size(),
                    this->core_material.size(),
                    this->shell_material.size(),
                    this->n_medium.size()
                };
            }
    };

    class Scatterer: public ScatteringProperties
    {
        public:
            double core_diameter;
            double shell_width;
            double shell_diameter;
            complex128 core_index;
            complex128 shell_index;
            double n_medium;
            double x_core;
            double x_shell;
            std::vector<complex128> an, bn;

            pybind11::array_t<complex128> get_an_py(){ return vector_to_numpy(this->an, {max_order}); }
            pybind11::array_t<complex128> get_bn_py(){ return vector_to_numpy(this->bn, {max_order}); }

            std::vector<complex128> get_an() const { return an; };
            std::vector<complex128> get_bn() const { return bn; };

            Scatterer(double wavelength, double amplitude, double core_diameter, double shell_width,
            complex128 core_index, complex128 shell_index, double n_medium, std::vector<complex128> jones, size_t max_order=0)
            : ScatteringProperties(wavelength, jones, amplitude), core_diameter(core_diameter), shell_width(shell_width),
            core_index(core_index), shell_index(shell_index), n_medium(n_medium)
            {
                this->shell_diameter = this->core_diameter + this->shell_width;
                initialize(max_order);
            }

            Scatterer(double core_diameter, double shell_width, complex128 core_index, complex128 shell_index,
                double n_medium, SOURCE::BaseSource &source, size_t max_order=0)
            : ScatteringProperties(source), core_diameter(core_diameter), shell_width(shell_width),
            core_index(core_index), shell_index(shell_index), n_medium(n_medium)
            {
                this->shell_diameter = this->core_diameter + this->shell_width;
                initialize(max_order);
            }

            void apply_medium(){
                this->core_index /= this->n_medium;
                this->shell_index /= this->n_medium;
                this->core_diameter *= this->n_medium;
                this->shell_width *= this->n_medium;
                this->shell_diameter *= this->n_medium;
            }

            void initialize(size_t &max_order);
            void compute_size_parameter();
            void compute_max_order(size_t &max_order);
            void compute_area();
            double get_g();
            double get_Qsca();
            double get_Qext();
            double get_Qback();

            void compute_an_bn();
            std::tuple<std::vector<complex128>, std::vector<complex128>> compute_s1s2(const std::vector<double> &phi);
    };
}


// -
