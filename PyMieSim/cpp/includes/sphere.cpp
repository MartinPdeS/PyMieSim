#ifndef SPHERE_H
#define SPHERE_H

  #include "base_scatterer.cpp"
  #include "numpy_interface.cpp"
  #include "VSH.cpp"

  namespace SPHERE
  {
    struct State
    {
      double diameter;
      complex128 index;
      double n_medium;
      State(){}
      State(double &diameter, complex128 &index, double &n_medium)
                  : diameter(diameter), index(index), n_medium(n_medium){}

      State(double &diameter, complex128 &&index, double &n_medium)
                  : diameter(diameter), index(index), n_medium(n_medium){}

      void apply_medium()
      {
        this->index    /= this->n_medium;
        this->diameter *= this->n_medium;
      }
    };


    class Scatterer: public ScatteringProperties
    {

        public:
            CVector an, bn, cn, dn;

            State state;

            Cndarray get_an_py() { return vector_to_ndarray(an, {max_order}); }
            Cndarray get_bn_py() { return vector_to_ndarray(bn, {max_order}); }
            Cndarray get_cn_py() { return vector_to_ndarray(cn, {max_order}); }
            Cndarray get_dn_py() { return vector_to_ndarray(dn, {max_order}); }

            CVector get_an() { return an; };
            CVector get_bn() { return bn; };
            CVector get_cn() { return cn; };
            CVector get_dn() { return dn; };

            Scatterer(){}
            Scatterer(double &wavelength,
                      double &amplitude,
                      double &diameter,
                      complex128 &index,
                      double &n_medium,
                      CVector &jones_vector,
                      size_t max_order = 0 )
                      : state(diameter, index, n_medium), ScatteringProperties(wavelength, jones_vector, amplitude)
                      {
                        initialize(max_order);
                      }

            Scatterer(State &state, SOURCE::State &source, size_t max_order = 0)
                      : state(state), ScatteringProperties(source)
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


        void compute_max_order(size_t &max_order)
        {
          if (max_order == 0)
            this->max_order = (size_t) (2 + this->size_parameter + 4 * pow(this->size_parameter, 1./3.)) + 16;
          else
            this->max_order = max_order;
        }

        void compute_size_parameter()
        {
          this->size_parameter = PI * state.diameter / source.wavelength;
        }

        void compute_area()
        {
          this->area = PI * pow(state.diameter/2.0, 2);
        }

        void compute_an_bn()
        {
          an = CVector(max_order);
          bn = CVector(max_order);

          complex128 x = size_parameter,
                     m = state.index,
                     mx = m * x,
                     _da, _db, _gsx, _gs1x, _px, _chx, _p1x, _ch1x, _p2x, _ch2x;

          size_t nmx = std::max( max_order, (size_t) std::abs(mx) ) + 16;
          CVector Dn  = VSH::SPHERICAL::compute_dn(nmx, mx);

          double n;

          _p1x  = sin(x);
          _ch1x = cos(x);

          for (size_t i = 1; i < max_order+1; ++i)
            {
                n = (double) i;
                _px =  x * compute_jn(n, x);
                _chx = -x * compute_yn(n, x);

                _p2x = _px;
                _ch2x = _chx;

                _gsx =  _px  - 1.*JJ * _chx;
                _gs1x =  _p1x - 1.*JJ * _ch1x;

                _da = Dn[i]/m + n/x;
                _db = Dn[i]*m + n/x;

                an[i-1] = (_da * _px - _p1x) / (_da * _gsx - _gs1x);
                bn[i-1] = (_db * _px - _p1x) / (_db * _gsx - _gs1x);

                _p1x  = _p2x;
                _ch1x = _ch2x;
            }
        }

      void compute_cn_dn()
      {
        cn = CVector(max_order);
        dn = CVector(max_order);

        complex128 x = size_parameter,
                   m = state.index,
                   z = m * x;

        size_t nmx = std::max( max_order, (size_t) std::abs(z) ) + 16;

        CVector Cnx = CVector(nmx),
                Cnn, jnx, jnmx, yx, hx, b1x, y1x, hn1x, ax, ahx, numerator,
                c_denominator, d_denominator;

        b1x.push_back( +sin(x) / x );
        y1x.push_back( -cos(x) / x );

        for (double i = nmx; i > 1; i--)
        {
          Cnx[i-2] = i - z*z/(Cnx[i-1] + i);
        }

        for (size_t i = 0; i < max_order; i++)
        {
          double n = (double) i;
          Cnn.push_back(Cnx[i]);
          jnx.push_back(compute_jn( n+1, x ));

          jnmx.push_back(1. / ( compute_jn(n+1, z )));
          yx.push_back(compute_yn(n+1, x ));
          hx.push_back(jnx[i] + JJ * yx[i]);

          b1x.push_back(jnx[i]);
          y1x.push_back(yx[i]);
          hn1x.push_back(b1x[i] + JJ * y1x[i]);

          ax.push_back(x * b1x[i] - ( n+1 ) * jnx[i]);
          ahx.push_back(x * hn1x[i] - ( n+1 ) * hx[i]);

          numerator.push_back( jnx[i] * ahx[i] - hx[i] * ax[i] );
          c_denominator.push_back( ahx[i] - hx[i] * Cnn[i] );
          d_denominator.push_back( m * m * ahx[i] - hx[i] * Cnn[i] );
          cn[i] = jnmx[i] * numerator[i] / c_denominator[i] ;
          dn[i] = jnmx[i] * m * numerator[i] / d_denominator[i] ;
        }
      }


      double get_g()
      {
          double value = 0;

          CVector an = get_an(), bn = get_bn();

          for(size_t it = 0; it < max_order-1; ++it)
          {
              double n = (double) it + 1;
              value += ( n * (n + 2.) / (n + 1.) ) * std::real(an[it] * std::conj(an[it+1]) + bn[it] * std::conj(bn[it+1]) );
              value += ( (2. * n + 1. ) / ( n * (n + 1.) ) )  * std::real( an[it] * std::conj(bn[it]) );
          }

          return value * 4. / ( get_Qsca() * pow(size_parameter, 2) );
      }

    std::tuple<CVector, CVector> compute_s1s2(const DVector &phi)
    {
      CVector an = get_an(), bn = get_bn();

      CVector S1(phi.size(), 0.0), S2(phi.size(), 0.0);

      DVector prefactor = get_prefactor();

      DVector Mu; Mu.reserve(phi.size());

      for (double phi : phi)
          Mu.push_back( cos( phi-PI / 2.0 ) );


      for (uint i = 0; i < phi.size(); i++)
      {
          auto [pin, taun] = VSH::SPHERICAL::MiePiTau( Mu[i], max_order);

          for (uint m = 0; m < max_order ; m++){
              S1[i] += prefactor[m] * ( an[m] * pin[m] +  bn[m] * taun[m] );
              S2[i] += prefactor[m] * ( an[m] * taun[m] + bn[m] * pin[m]  );
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

        value = pow( std::abs(value), 2. ) / pow( size_parameter, 2. );
        return std::abs(value);
    }
  };
}


#endif
// -
