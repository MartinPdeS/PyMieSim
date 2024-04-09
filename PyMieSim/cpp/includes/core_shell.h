#pragma once

#include "base_scatterer.cpp"
#include "numpy_interface.cpp"
#include "VSH.cpp"

  namespace CORESHELL
  {

    struct State
    {
      double core_diameter, shell_width, shell_diameter;
      complex128 core_index, shell_index;
      double n_medium, x_core, x_shell;

      State() = default;

        State(double core_diameter, double shell_width, complex128 core_index,
        complex128 shell_index, double n_medium)
        : core_diameter(core_diameter), shell_width(shell_width), core_index(core_index),
        shell_index(shell_index), n_medium(n_medium) {
            this->shell_diameter = this->core_diameter + this->shell_width;
        }

        void apply_medium(){
            this->core_index /= this->n_medium;
            this->shell_index /= this->n_medium;
            this->core_diameter *= this->n_medium;
            this->shell_width *= this->n_medium;
            this->shell_diameter *= this->n_medium;
        }
    };

    class Scatterer: public ScatteringProperties
    {
        public:
            CVector an, bn;

            State state;

            Cndarray get_an_py(){ return vector_to_ndarray(this->an, {max_order}); }
            Cndarray get_bn_py(){ return vector_to_ndarray(this->bn, {max_order}); }

            std::vector<complex128> get_an(){ return an; };
            std::vector<complex128> get_bn(){ return bn; };

            Scatterer(double wavelength, double amplitude, double core_diameter, double shell_width,
            complex128 core_index, complex128 shell_index, double n_medium, std::vector<complex128> jones, size_t max_order=0)
            : ScatteringProperties(wavelength, jones, amplitude), state(core_diameter, shell_width, core_index, shell_index, n_medium){
                initialize(max_order);
            }

            Scatterer(State &state, SOURCE::State &source, size_t max_order=0)
            : ScatteringProperties(source), state(state){
                initialize(max_order);
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
            std::tuple<CVector, CVector> compute_s1s2(const DVector &Phi);
    };
}


// -
