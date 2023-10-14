#ifndef CORESHELL_H
#define CORESHELL_H

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
      State(){}
      State(double &core_diameter,
            double &shell_width,
            complex128 &core_index,
            complex128& shell_index,
            double &n_medium
            ) : core_diameter(core_diameter),
                shell_width(shell_width),
                core_index(core_index),
                shell_index(shell_index),
                n_medium(n_medium)
      {
        this->shell_diameter = this->core_diameter + this->shell_width;
      }

      void apply_medium()
      {
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

      Scatterer(double &Wavelength,
                double &Amplitude,
                double &core_diameter,
                double &shell_width,
                complex128 &core_index,
                complex128 &shell_index,
                double &n_medium,
                std::vector<complex128> &Jones,
                size_t max_order=0 )
                : ScatteringProperties(Wavelength, Jones, Amplitude), state(core_diameter, shell_width, core_index, shell_index, n_medium)
                {
                  initialize(max_order);
                }

      Scatterer(State &state, SOURCE::State &source, size_t max_order=0) : state(state), ScatteringProperties(source)
      {
        initialize(max_order);
      }

      void initialize(size_t &max_order)
      {
        state.apply_medium();
        compute_size_parameter();
        compute_max_order(max_order);
        compute_area();
        compute_an_bn();
      }

      void compute_size_parameter()
      {
        size_parameter = source.k * state.shell_diameter / 2;
        state.x_shell = source.k * state.shell_diameter / 2.0;
        state.x_core = source.k * state.core_diameter / 2.0;
      }

      void compute_max_order(size_t &max_order)
      {
        if (max_order == 0)
          this->max_order = (size_t) (2 + size_parameter + 4 * pow(size_parameter,1./3.));
        else
          this->max_order = max_order;
      }

      void compute_area()
      {
         this->area = PI * pow(state.shell_diameter/2.0, 2);
      }

      double get_g()
      {
          double Qsca = get_Qsca(),
                 value = 0;

          for(size_t it = 0; it < max_order-1; ++it)
          {
              double n = (double) it + 1;
              value += ( n * (n + 2.) / (n + 1.) ) * std::real(an[it] * std::conj(an[it+1]) + bn[it] * std::conj(bn[it+1]) );
              value += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( an[it] * std::conj(bn[it]) );
          }

          return value * 4. / ( Qsca * pow(size_parameter, 2) );
      }

    void compute_an_bn()
    {
      an = CVector(max_order);
      bn = CVector(max_order);

      complex128 m = state.shell_index / state.core_index,
                 u = state.core_index * state.x_core,
                 v = state.shell_index * state.x_core,
                 w = state.shell_index * state.x_shell;

      complex128 sv = sqrt(0.5 * PI * v),
                 sw = sqrt(0.5 * PI * w),
                 sy = sqrt(0.5 * PI * state.x_shell);

      size_t mx   = (size_t) std::max( abs( state.core_index * state.x_shell ), abs( state.shell_index*state.x_shell ) ),
             nmx  = (size_t) ( std::max( max_order, mx ) + 16. )  ;

      CVector pv, pw, py, chv, chw, chy, p1y, ch1y, gsy, gs1y;

      p1y. push_back( sin( state.x_shell ) ) ;
      ch1y.push_back( cos( state.x_shell ) ) ;

      for (size_t i=0; i<max_order+1; i++)
      {
        double nu = i + 1.5 ;
        pw.push_back( sw * compute_Jn(nu, w) );
        pv.push_back( sv * compute_Jn(nu, v) );
        py.push_back( sy * compute_Jn(nu, state.x_shell) );

        chv.push_back( -sv * compute_Yn(nu, v) );
        chw.push_back( -sw * compute_Yn(nu, w) );
        chy.push_back( -sy * compute_Yn(nu, state.x_shell) );

        p1y.push_back ( py[i]  );
        ch1y.push_back( chy[i] );
        gsy.push_back ( py[i]  - JJ * chy[i]  );
        gs1y.push_back ( p1y[i] - JJ * ch1y[i] );
      }

      CVector Du = CVector(nmx, 0.),
              Dv = CVector(nmx, 0.),
              Dw = CVector(nmx, 0.);

      for (int i = nmx-1; i > 1; i--)
      {
        Du[i-1] = (double)i / u -1. / (Du[i] + (double)i / u);
        Dv[i-1] = (double)i / v -1. / (Dv[i] + (double)i / v);
        Dw[i-1] = (double)i / w -1. / (Dw[i] + (double)i / w);
      }

      Du.erase(Du.begin());
      Dv.erase(Dv.begin());
      Dw.erase(Dw.begin());

      CVector uu, vv, fv, dns, gns, a1, b1;
      for (size_t i=0; i<max_order; i++)
      {
        double n = (double) (i+1);
        uu.push_back ( m * Du[i] - Dv[i]  );
        vv.push_back ( Du[i] / m - Dv[i] );
        fv.push_back ( pv[i] / chv[i]    );
        dns.push_back( ( ( uu[i] * fv[i] / pw[i] ) / ( uu[i] * ( pw[i] - chw[i] * fv[i] ) + ( pw[i] / pv[i] ) / chv[i] ) ) + Dw[i] );
        gns.push_back( ( ( vv[i] * fv[i] / pw[i] ) / ( vv[i] * ( pw[i] - chw[i] * fv[i] ) + ( pw[i] / pv[i] ) / chv[i] ) ) + Dw[i] );
        a1.push_back ( dns[i] / state.shell_index + n / state.x_shell );
        b1.push_back ( state.shell_index * gns[i] + n / state.x_shell );
        an[i] = ( py[i] * a1[i] - p1y[i] ) / ( gsy[i] * a1[i] - gs1y[i] ) ;
        bn[i] = ( py[i] * b1[i] - p1y[i] ) / ( gsy[i] * b1[i] - gs1y[i] ) ;
      }
    }


    std::tuple<CVector, CVector> compute_s1s2(const DVector &Phi)
    {
      CVector S1(Phi.size(), 0.0), S2(Phi.size(), 0.0);

      DVector Prefactor = get_prefactor();

      DVector Mu; Mu.reserve(Phi.size());

      for (double phi : Phi)
          Mu.push_back( cos( phi-PI/2.0 ) );


      for (uint i = 0; i < Phi.size(); i++)
      {
          auto [pin, taun] = VSH::SPHERICAL::MiePiTau( Mu[i], max_order);

          for (uint m = 0; m < max_order ; m++){
              S1[i]    += Prefactor[m] * ( an[m] * pin[m] +  bn[m] * taun[m] );
              S2[i]    += Prefactor[m] * ( an[m] * taun[m] + bn[m] * pin[m]  );
            }
      }

      return std::make_tuple(S1, S2)  ;
    }


    double get_Qsca()
    {
        CVector an = get_an(), bn = get_bn();

        double value = 0;
        for(size_t it = 0; it < max_order; ++it)
        {
             double n = (double) it + 1;
             value += (2.* n + 1.) * ( pow( std::abs(an[it]), 2) + pow( std::abs(bn[it]), 2)  );

        }
        return value * 2. / pow( size_parameter, 2.);
    }


    double get_Qext()
    {
      CVector an = get_an(), bn = get_bn();

      double value = 0;
      for(size_t it = 0; it < max_order; ++it)
      {
           double n = (double) it + 1;
           value += (2.* n + 1.) * std::real( an[it] + bn[it] );

      }
      return value * 2. / pow( size_parameter, 2.);
    }


    double get_Qback()
    {
        CVector an = get_an(), bn = get_bn();

        complex128 value = 0;

        for(size_t it = 0; it < max_order-1; ++it)
        {
          double n = (double) it + 1;
          value += (2. * n + 1) * pow(-1., n) * ( an[it] - bn[it] ) ;
        }

        value *= pow( std::abs(value), 2. ) / pow( size_parameter, 2. );
        return std::abs(value);
    }
  };
}


#endif


// -
