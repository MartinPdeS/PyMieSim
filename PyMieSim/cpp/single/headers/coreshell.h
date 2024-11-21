#pragma once

#include "single/includes/base_spherical_scatterer.cpp"
#include "utils/VSH.cpp"


  namespace CORESHELL
  {
    class Scatterer: public BaseSphericalScatterer
    {
        public:
            double core_diameter;
            double shell_thickness;
            double shell_diameter;
            complex128 core_index;
            complex128 shell_index;
            double x_core;
            double x_shell;
            std::vector<size_t> indices;

            Scatterer(double core_diameter, double shell_thickness, complex128 core_index, complex128 shell_index,
                double medium_index, SOURCE::BaseSource &source, size_t max_order = 0)
            : BaseSphericalScatterer(source, max_order, medium_index), core_diameter(core_diameter), shell_thickness(shell_thickness),
            core_index(core_index), shell_index(shell_index)
            {
                this->shell_diameter = this->core_diameter + this->shell_thickness;
                this->compute_area();
                this->compute_size_parameter();
                this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
                this->apply_medium();
                this->compute_an_bn();
            }

            void compute_size_parameter() override {
                this->size_parameter = source.wavenumber * this->shell_diameter / 2 * this->medium_index;
                this->size_parameter_squared = pow(this->size_parameter, 2);
                this->x_shell = source.wavenumber * this->shell_diameter / 2.0 * this->medium_index;
                this->x_core = source.wavenumber * this->core_diameter / 2.0 * this->medium_index;
            }
            void compute_area() override {

                this->area = PI * pow(this->shell_diameter/2.0, 2);
            }

            void apply_medium();
            void compute_cn_dn();
            void compute_an_bn();
            void compute_max_order(size_t max_order);
    };
}


// -
