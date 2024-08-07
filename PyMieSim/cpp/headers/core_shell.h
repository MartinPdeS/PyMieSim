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

            Scatterer(double core_diameter, double shell_width, complex128 core_index, complex128 shell_index,
                double medium_index, SOURCE::BaseSource &source, size_t max_order = 0)
            : BaseSphericalScatterer(source, max_order, medium_index), core_diameter(core_diameter), shell_width(shell_width),
            core_index(core_index), shell_index(shell_index)
            {
                this->shell_diameter = this->core_diameter + this->shell_width;
                this->compute_area();
                this->compute_size_parameter();
                this->max_order = (max_order == 0) ? this->get_wiscombe_criterion(this->size_parameter) : max_order;
                this->apply_medium();
                this->compute_an_bn();
            }

            void compute_size_parameter() override {
                size_parameter = source.k * this->shell_diameter / 2;
                this->x_shell = source.k * this->shell_diameter / 2.0;
                this->x_core = source.k * this->core_diameter / 2.0;
            }
            void compute_area() override {

                this->area = PI * pow(this->shell_diameter/2.0, 2);
            }

            void apply_medium();
            void compute_cn_dn();
            void compute_an_bn();
            void compute_max_order(size_t max_order);
    };

    class Set {
        public:
            std::vector<double> core_diameter;
            std::vector<double> shell_width;
            std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> core;
            std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> shell;
            std::variant<std::vector<double>, std::vector<std::vector<double>>> medium;
            std::vector<size_t> shape;

            Set() = default;

            Set(
                const std::vector<double>& core_diameter,
                const std::vector<double>& shell_width,
                std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> core_param,
                std::variant<std::vector<complex128>, std::vector<std::vector<complex128>>> shell_param,
                std::variant<std::vector<double>, std::vector<std::vector<double>>> medium_param
            ) : core_diameter(core_diameter), shell_width(shell_width),
                core(core_param), shell(shell_param), medium(medium_param) {
                update_shape();
            }

            void update_shape() {
                shape.clear();
                shape.push_back(core_diameter.size());
                shape.push_back(shell_width.size());

                // Safely obtaining sizes
                shape.push_back(
                    std::holds_alternative<std::vector<complex128>>(core) ?
                    std::get<std::vector<complex128>>(core).size() :
                    std::get<std::vector<std::vector<complex128>>>(core).size()
                );

                shape.push_back(
                    std::holds_alternative<std::vector<complex128>>(shell) ?
                    std::get<std::vector<complex128>>(shell).size() :
                    std::get<std::vector<std::vector<complex128>>>(shell).size()
                );

                shape.push_back(
                    std::holds_alternative<std::vector<double>>(medium) ?
                    std::get<std::vector<double>>(medium).size() :
                    std::get<std::vector<std::vector<double>>>(medium).size()
                );
            }

            Scatterer to_object(size_t wl, size_t cd, size_t sw, size_t ci, size_t si, size_t mi, SOURCE::BaseSource &source) const {
                complex128 core_value;
                complex128 shell_value;
                double medium_value;

                if (std::holds_alternative<std::vector<std::vector<complex128>>>(core)) {
                    const auto& materials = std::get<std::vector<std::vector<complex128>>>(core);
                    core_value = materials[ci][wl];
                } else {
                    const auto& indices = std::get<std::vector<complex128>>(core);
                    core_value = indices[ci];
                }

                if (std::holds_alternative<std::vector<std::vector<complex128>>>(shell)) {
                    const auto& materials = std::get<std::vector<std::vector<complex128>>>(shell);
                    shell_value = materials[si][wl];
                } else {
                    const auto& indices = std::get<std::vector<complex128>>(shell);
                    shell_value = indices[si];
                }

                if (std::holds_alternative<std::vector<std::vector<double>>>(medium)) {
                    const auto& mat = std::get<std::vector<std::vector<double>>>(medium);
                    medium_value = mat[mi][wl];
                } else {
                    const auto& indices = std::get<std::vector<double>>(medium);
                    medium_value = indices[mi];
                }

                return Scatterer(core_diameter[cd], shell_width[sw], core_value, shell_value, medium_value, source);
            }

    };

}


// -
