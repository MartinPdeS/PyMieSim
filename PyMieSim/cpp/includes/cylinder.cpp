#ifndef CYLINDER_H
#define CYLINDER_H

#include "base_scatterer.cpp"

namespace CYLINDER
{

  struct State
  {
    double
      diameter;

    complex128
      index;

    double
      n_medium;

    State(){}
    State(
      const double &diameter,
      const complex128 &index,
      const double &n_medium
    )
      : diameter(diameter),
        index(index),
        n_medium(n_medium){}

    void apply_medium()
    {
      this->index /= this->n_medium;
      this->diameter *= this->n_medium;
    }
  };


  class Scatterer: public ScatteringProperties
  {
      public:

          std::vector<complex128>
            a1n,
            b1n,
            a2n,
            b2n;

          State state;

          Cndarray get_a1n_py() { return vector_to_ndarray(a1n, {max_order}); }
          Cndarray get_b1n_py() { return vector_to_ndarray(b1n, {max_order}); }
          Cndarray get_a2n_py() { return vector_to_ndarray(a2n, {max_order}); }
          Cndarray get_b2n_py() { return vector_to_ndarray(b2n, {max_order}); }

          std::vector<complex128> get_a1n() { return a1n; };
          std::vector<complex128> get_b1n() { return b1n; };
          std::vector<complex128> get_a2n() { return a2n; };
          std::vector<complex128> get_b2n() { return b2n; };

          double get_diameter(){return state.diameter;}
          complex128 get_index(){return state.index;}
          double get_wavelength(){return source.wavelength;}
          double get_k(){return source.k;}
          std::vector<complex128> get_jones_vector(){return source.jones_vector;}

          Scatterer(
            double wavelength,
            double amplitude,
            double diameter,
            complex128 index,
            double n_medium,
            std::vector<complex128> jones_vector,
            size_t max_order = 0
          )
            : ScatteringProperties(wavelength, jones_vector, amplitude),
              state(diameter, index, n_medium)
            {
              initialize(max_order);
            }

          Scatterer(
            State &state,
            SOURCE::State &source,
            size_t max_order = 0)
            : ScatteringProperties(source),
              state(state)
            {
              initialize(max_order);
            }

      void initialize(size_t &max_order)
      {
        // state.apply_medium();
        this->compute_size_parameter();
        this->compute_max_order(max_order);
        this->compute_area();
        this->compute_an_bn();

      }

      void compute_max_order(size_t &max_order)
      {
        if (max_order == 0)
          this->max_order  = (size_t) (2 + this->size_parameter + 4 * pow(this->size_parameter, 1./3.)) + 16;
        else
          this->max_order = max_order;
      }

      void compute_size_parameter()
      {
        this->size_parameter = PI * state.diameter / source.wavelength;
      }

      void compute_area()
      {
        this->area = state.diameter;
      }

      double get_g()
      {
          return get_g_with_fields(1000, 1.0);
      }

      double get_Qsca()
      {
          std::vector<complex128>
            a1n = get_a1n(),
            b1n = get_b1n(),
            a2n = get_a2n(),
            b2n = get_b2n();

          complex128 Qsca1=0, Qsca2=0;

          for(size_t it = 1; it < max_order; it++)
          {
            Qsca1 +=  pow( std::abs(a1n[it]), 2) + pow( std::abs(b1n[it]), 2) ;
            Qsca2 +=  pow( std::abs(a2n[it]), 2) + pow( std::abs(b2n[it]), 2) ;
          }

          Qsca1 =  2. / size_parameter * ( 2.0 * Qsca1 + pow( abs(b1n[0]), 2 ) );
          Qsca2 =  2. / size_parameter * ( 2.0 * Qsca2 + pow( abs(a2n[0]), 2 ) );

          return process_polarization(Qsca1, Qsca2);
      }

      double get_Qext()
      {
        std::vector<complex128>
          a1n = get_a1n(),
          b1n = get_b1n(),
          a2n = get_a2n(),
          b2n = get_b2n();

        complex128
          Qext1 = 0,
          Qext2 = 0;

        for(size_t it = 1; it < max_order; ++it)
        {
          Qext1 += b1n[it];
          Qext2 += a2n[it];
        }

        Qext1 = 2. / size_parameter * std::real( b1n[0] + 2.0 * Qext1 );
        Qext2 = 2. / size_parameter * std::real( a1n[0] + 2.0 * Qext2 );

        return this->process_polarization(Qext1, Qext2);
      }

      double process_polarization(complex128 &Value0, complex128& value1)
      {
        if (source.is_polarized == false)
            return 0.5 * ( abs( value1 ) + abs( Value0 ) );
        else
            return abs( value1 ) * pow(abs(source.jones_vector[0]), 2) + abs( Value0 ) * pow(abs(source.jones_vector[1]), 2);
      }

      void compute_an_bn()
      {
        a1n = std::vector<complex128>(max_order);
        b1n = std::vector<complex128>(max_order);

        a2n = std::vector<complex128>(max_order);
        b2n = std::vector<complex128>(max_order);

        double
          x = size_parameter;

        complex128
          m = state.index / state.n_medium,
          z = m * x,
          numerator,
          denominator;

        std::vector<complex128>
          J_z(max_order+1),
          J_z_p(max_order+1),
          J_x(max_order+1),
          J_x_p(max_order+1),
          H_x(max_order+1),
          H_x_p(max_order+1);

        for (size_t n = 0; n < max_order+1; ++n)
        {
          double
            nd = (double) n ;

          J_z[n] = compute_Jn(nd, z);
          J_z_p[n] = compute_Jn_p(nd, z);
          J_x[n] = compute_Jn(nd, x);
          J_x_p[n] = compute_Jn_p(nd, x);
          H_x[n] = compute_H1(nd, x);
          H_x_p[n] = compute_H1_p(nd, x);
        }

        for (size_t n = 0; n < (size_t) max_order; n++)
        {
          numerator = m * J_z[n] * J_x_p[n] - J_z_p[n] * J_x[n];
          denominator = m * J_z[n] * H_x_p[n] - J_z_p[n] * H_x[n];
          a1n[n] = 0.0 ;
          a2n[n] = numerator/denominator ;

          numerator = J_z[n] * J_x_p[n] - m * J_z_p[n] * J_x[n];
          denominator = J_z[n] * H_x_p[n] - m * J_z_p[n] * H_x[n];
          b1n[n] = numerator/denominator ;
          b2n[n] = 0.0 ;
        }
      }

      std::tuple<std::vector<complex128>, std::vector<complex128>>
      compute_s1s2(const std::vector<double> &phi)
      {
        std::vector<complex128>
          T1(phi.size()),
          T2(phi.size());

        for (uint i = 0; i < phi.size(); i++){
            T1[i] = b1n[0];
            T2[i] = a2n[0];
            for (size_t order = 1; order < max_order ; order++){
                T1[i] += 2.0 * b1n[order] * cos(order * (PI - (phi[i] + PI/2.0) ) );
                T2[i] += 2.0 * a2n[order] * cos(order * (PI - (phi[i] + PI/2.0) ) );
              }
        }

        return std::make_tuple( T1, T2 )  ;
      }

  };

}

#endif
