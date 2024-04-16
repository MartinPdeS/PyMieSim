#pragma once

#include "base_spherical_scatterer.cpp"
#include "numpy_interface.cpp"
#include "VSH.cpp"

  namespace CORESHELL
  {
    class Scatterer: public BaseSphericalScatterer
    {
        public:
            double core_diameter;
            double shell_width;
            double shell_diameter;
            complex128 core_index;
            complex128 shell_index;
            double x_core;
            double x_shell;

            Scatterer(double wavelength, double amplitude, double core_diameter, double shell_width,
            complex128 core_index, complex128 shell_index, double n_medium, std::vector<complex128> jones, size_t max_order)
            : BaseSphericalScatterer(wavelength, jones, amplitude, n_medium), core_diameter(core_diameter), shell_width(shell_width),
            core_index(core_index), shell_index(shell_index)
            {
                this->shell_diameter = this->core_diameter + this->shell_width;
                apply_medium();
                compute_size_parameter();
                this->max_order = max_order;
                compute_area();
                compute_an_bn();
            }

            Scatterer(double wavelength, double amplitude, double core_diameter, double shell_width,
            complex128 core_index, complex128 shell_index, double n_medium, std::vector<complex128> jones)
            : BaseSphericalScatterer(wavelength, jones, amplitude, n_medium), core_diameter(core_diameter), shell_width(shell_width),
            core_index(core_index), shell_index(shell_index)
            {
                this->shell_diameter = this->core_diameter + this->shell_width;
                apply_medium();
                compute_size_parameter();
                this->max_order = get_wiscombe_criterion(this->size_parameter);
                compute_area();
                compute_an_bn();
            }

            Scatterer(double core_diameter, double shell_width, complex128 core_index, complex128 shell_index,
                double n_medium, SOURCE::BaseSource &source, size_t max_order)
            : BaseSphericalScatterer(source, n_medium), core_diameter(core_diameter), shell_width(shell_width),
            core_index(core_index), shell_index(shell_index)
            {
                this->shell_diameter = this->core_diameter + this->shell_width;
                apply_medium();
                compute_size_parameter();
                this->max_order = max_order;
                compute_area();
                compute_an_bn();
            }

            Scatterer(double core_diameter, double shell_width, complex128 core_index, complex128 shell_index,
                double n_medium, SOURCE::BaseSource &source)
            : BaseSphericalScatterer(source, n_medium), core_diameter(core_diameter), shell_width(shell_width),
            core_index(core_index), shell_index(shell_index)
            {
                this->shell_diameter = this->core_diameter + this->shell_width;
                apply_medium();
                compute_size_parameter();
                this->max_order = get_wiscombe_criterion(this->size_parameter);
                compute_area();
                compute_an_bn();
            }



            void apply_medium();
            void compute_cn_dn();
            void compute_an_bn();
            void compute_max_order(size_t max_order);
            void compute_size_parameter();
            void compute_area();
    };

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

            Scatterer to_object(size_t wl, size_t cd, size_t sw, size_t ci, size_t si, size_t mi, SOURCE::BaseSource &source) const
            {
                if (core_is_material && shell_is_material)
                    return Scatterer(
                        this->core_diameter[cd],
                        this->shell_width[sw],
                        this->core_material[ci][wl],
                        this->shell_material[si][wl],
                        this->n_medium[mi],
                        source
                    );

                if (core_is_material && !shell_is_material)
                    return Scatterer(
                        this->core_diameter[cd],
                        this->shell_width[sw],
                        this->core_material[ci][wl],
                        this->shell_index[si],
                        this->n_medium[mi],
                        source
                    );

                if (!core_is_material && shell_is_material)
                    return Scatterer(
                        this->core_diameter[cd],
                        this->shell_width[sw],
                        this->core_index[ci],
                        this->shell_material[si][wl],
                        this->n_medium[mi],
                        source
                    );

                if (!core_is_material && !shell_is_material)
                    return Scatterer(
                        this->core_diameter[cd],
                        this->shell_width[sw],
                        this->core_index[ci],
                        this->shell_index[si],
                        this->n_medium[mi],
                        source
                    );






            }
    };
}


// -
