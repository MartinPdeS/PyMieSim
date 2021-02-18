
#include <vector>
#include <complex>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>
typedef std::complex<double> complex128;
typedef std::vector<complex128> iVec;
typedef std::vector<double> Vec;
#define PI 3.14159265


std::tuple<Vec, Vec, Vec>
Cart2sph(Vec X, Vec Y, Vec Z){
  Vec R, PHI, THETA;
  for (long unsigned int i = 0; i < X.size(); i++){
    R.push_back( sqrt( pow( X[i],2) + pow(Y[i],2) + pow(Z[i],2) ) );
    PHI.push_back( atan( Y[i]/X[i] ) );
    THETA.push_back( atan( sqrt( pow(Y[i],2) / pow(X[i],2) ) ) );
  }
  return std::make_tuple(R, PHI, THETA);
}


std::tuple<Vec, Vec, Vec>
Fibonacci_sphere(int samples,
                 double maxAngle){


    Vec X, Y , Z;

    Vec R, THETA, PHI;
    double phi = PI * (3. - sqrt(5.));


    double SolidAngle = abs( 2*PI * (cos(maxAngle) - 1));

    double ratio = 4*PI / SolidAngle;

    int TrueSampling = (int)( samples * ratio );

    for (long unsigned int i = 0; i < TrueSampling; i++){

        double y = 1 - (i / (double)(TrueSampling - 1)) * 2 ;
        double radius = sqrt(1 - y * y);

        double theta = phi * i;

        double x = cos(theta) * radius;
        double z = sin(theta) * radius;

        X.push_back(x);
        Y.push_back(1.0);
        Z.push_back(1.0);

        if (i >= samples - 1) break;
      }

    std::tie(R, PHI, THETA) = Cart2sph(X,Y,Z);

    X.clear();
    Y.clear();
    Z.clear();

    X.shrink_to_fit();
    Y.shrink_to_fit();
    Z.shrink_to_fit();

    return std::make_tuple(R, PHI, THETA);


}




double
Legendre(std::size_t n, double x)
{
  return boost::math::legendre_p(n, 1, x);
}

double
Legendre_p(std::size_t n, double x)
{
  return -boost::math::legendre_p(n, 1, x)/(std::sqrt(1.0 - x*x));
}

double
Legendre_pp(std::size_t n, double x)
{
  return boost::math::legendre_p(n, 2, x)/(1.0 - x*x);
}

double
Pin(std::size_t n, double theta)
{
  return -1 * Legendre_p(n, cos(theta));
}


double
Taun(std::size_t n, double theta)
{
  return - cos(theta) * Legendre_p(n, cos(theta)) + sin(theta)*sin(theta) * Legendre_pp(n, cos(theta));
}



template <class T>
std::vector<T>
Riccati1(std::size_t n, T x)
{
  std::vector<T> temp;
  T b = boost::math::sph_bessel(n, x);
  T b_p = boost::math::sph_bessel_prime(n, x);
  temp.push_back(x * b);
  temp.push_back(b + x * b_p);

  return temp;
}


template <class T>
std::vector<T>
Riccati2(std::size_t n, T x)
{
  std::vector<T> temp;
  T b = boost::math::sph_neumann(n, x);
  T b_p = boost::math::sph_neumann_prime(n, x);
  temp.push_back(x * b);
  temp.push_back(b + x * b_p);

  return temp;
}


iVec
Riccati3(std::size_t n, double x)
{
  complex128 j (0., 1.0);
  iVec ret;
  Vec temp0, temp1;
  temp0 =  Riccati1(n,x);
  temp1 = Riccati2(n,x);


  ret.reserve(temp0.size());

  std::transform( begin(temp0),
                  end(temp0),
                  begin(temp1),
                  std::back_inserter(ret),
                  [](double r, double i) { return std::complex<double>(r, i); }
                 );

  return ret;

}


template <class T>
std::vector<T>
VecHankel(std::size_t n, T x)
{
  complex128 j (0., 1.0);
  iVec ret;
  Vec temp0, temp1;
  temp0 = boost::math::sph_bessel(n, x);
  temp1 = boost::math::sph_neumann(n, x);


  ret.reserve(temp0.size());

  std::transform( begin(temp0),
                  end(temp0),
                  begin(temp1),
                  std::back_inserter(ret),
                  [](double r, double i) { return std::complex<double>(r, i); }
                 );

  return ret;
}


template <class T>
std::vector<T>
VecHankel_p(std::size_t n, T x)
{
  complex128 j (0., 1.0);
  iVec ret;
  Vec temp0, temp1;
  temp0 = boost::math::sph_bessel_prime(n, x);
  temp1 = boost::math::sph_neumann_prime(n, x);


  ret.reserve(temp0.size());

  std::transform( begin(temp0),
                  end(temp0),
                  begin(temp1),
                  std::back_inserter(ret),
                  [](double r, double i) { return std::complex<double>(r, i); }
                 );

  return ret;
}





complex128
Hankel(std::size_t n, double x)
{
  complex128 j (0., 1.0);
  complex128 ret;
  complex128 temp0, temp1;
  temp0 = boost::math::sph_bessel(n, x);
  temp1 = boost::math::sph_neumann(n, x);

  return temp0 + j * temp1;

}


complex128
Hankel_p(std::size_t n, double x)
{
  complex128 j (0., 1.0);
  complex128 ret;
  complex128 temp0, temp1;
  temp0 = boost::math::sph_bessel_prime(n, x);
  temp1 = boost::math::sph_neumann_prime(n, x);

  return temp0 + j * temp1;

}



class VectorSphericalHarmonics1{
  public:

    int n;

    double k;

    VectorSphericalHarmonics1(int n, double Wavelength){
      this->n = n;
      this->k =2 * PI /Wavelength;
    }

    double Pin(std::size_t n, double x){return Pin(n, x);}

    double Taun(std::size_t n, double x){return Taun(n, x);}

    double ZFunc(std::size_t n, double x){ return boost::math::sph_bessel(n, x);}

    double ZFunc_p(std::size_t n, double x){ return boost::math::sph_bessel_prime(n, x);}

        std::tuple<iVec*, iVec*, iVec*>
        M_o1n(Vec r, Vec theta, Vec phi){
          iVec *ThetaComp = new iVec(r.size());
          iVec *PhiComp   = new iVec(r.size());
          iVec *RComp     = new iVec(r.size());

          for (long unsigned int i = 0; i < r.size(); i++){
            double p = this->k*r[i];
            ThetaComp->push_back(cos(phi[i]) * this->Pin(this->n, theta[i]) * this->ZFunc(this->n,p)  );
            PhiComp->push_back( -sin(phi[i]) * this->Taun(this->n, theta[i]) * this->ZFunc(this->n, p) );
            RComp->push_back(0.);
          }

          return std::make_tuple(RComp, PhiComp, ThetaComp);
        };

        std::tuple<iVec*, iVec*, iVec*>
        M_e1n(Vec r, Vec theta, Vec phi){
          iVec *ThetaComp = new iVec(r.size());
          iVec *PhiComp   = new iVec(r.size());
          iVec *RComp     = new iVec(r.size());

          for (long unsigned int i = 0; i < r.size(); i++){
            double p = this->k*r[i];
            ThetaComp->push_back(-sin(phi[i]) * this->Pin(this->n, theta[i]) * this->ZFunc(this->n, p) );
            PhiComp->push_back(  -cos(phi[i]) * this->Taun(this->n, theta[i]) * this->ZFunc(this->n, p) );
            RComp->push_back(0.);
          }

          return std::make_tuple(RComp, PhiComp, ThetaComp);
        };

        std::tuple<iVec*, iVec*, iVec*>
        N_o1n(Vec r, Vec theta, Vec phi){
          iVec *ThetaComp = new iVec(r.size());
          iVec *PhiComp   = new iVec(r.size());
          iVec *RComp     = new iVec(r.size());

          for (long unsigned int i = 0; i < r.size(); i++){
            double p = this->k*r[i];
            ThetaComp->push_back(sin(phi[i]) * this->Taun(this->n, theta[i]) * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p)/p )  );
            PhiComp->push_back(  cos(phi[i]) * this->Pin(this->n, theta[i])  * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p)/p )  );
            RComp->push_back(    sin(phi[i]) * this->n * (this->n+1) * sin(theta[i]) * Pin(this->n, theta[i]) * ZFunc(this->n, p)/p  );
          }

          return std::make_tuple(RComp, PhiComp, ThetaComp);
        };

        std::tuple<iVec*, iVec*, iVec*>
        N_e1n(Vec r, Vec theta, Vec phi){

          iVec *ThetaComp = new iVec(r.size());
          iVec *PhiComp   = new iVec(r.size());
          iVec *RComp     = new iVec(r.size());

          for (long unsigned int i = 0; i < r.size(); i++){

            double p = this->k*r[i];
            ThetaComp->push_back(cos(phi[i]) * this->Taun(this->n, theta[i]) * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p)/p )  );
            PhiComp->push_back(  sin(phi[i]) * this->Pin(this->n, theta[i])  * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p)/p )  );
            RComp->push_back(    cos(phi[i]) * this->n * (n+1) * sin(theta[i]) * Pin(this->n, theta[i]) * ZFunc(this->n, p)/p  );
          }

          return std::make_tuple(RComp, PhiComp, ThetaComp);
        };

    };



std::tuple<iVec*, iVec*, iVec*>
N3e1n(int n, double k, Vec r, Vec theta, Vec phi){

  iVec *ThetaComp = new iVec(r.size());
  iVec *PhiComp   = new iVec(r.size());
  iVec *RComp     = new iVec(r.size());

  for (long unsigned int i = 0; i < r.size(); i++){

    double p = k*r[i];

    (*ThetaComp)[i] =  cos(phi[i]) * Taun(n, theta[i]) * ( Hankel(n, p) + p * Hankel_p(n, p) ) / p ;

    (*PhiComp)[i]   = -sin(phi[i]) * Pin(n, theta[i])  * ( Hankel(n, p) + p * Hankel_p(n, p) ) / p ;

    (*RComp)[i]     =  cos(phi[i]) * (complex128)( n * (n+1) ) * sin(theta[i]) * Pin(n, theta[i]) * Hankel(n, p) / p ;

  }

  return std::make_tuple(RComp, PhiComp, ThetaComp);
}




std::tuple<iVec*, iVec*, iVec*>
M3o1n(int n, double k, Vec r, Vec theta, Vec phi){

  iVec *ThetaComp = new iVec(r.size());
  iVec *PhiComp   = new iVec(r.size());
  iVec *RComp     = new iVec(r.size());

  for (long unsigned int i = 0; i < r.size(); i++){

    double p = k*r[i];

    (*ThetaComp)[i] = ( cos(phi[i]) * Pin(n, theta[i]) * Hankel(n,p)  );

    (*PhiComp)[i]   = ( -sin(phi[i]) * Taun(n, theta[i]) * Hankel(n, p) );

    (*RComp)[i]     = (0.);

  }

  return std::make_tuple(RComp, PhiComp, ThetaComp);
}



std::tuple<iVec*, iVec*, iVec*>
N3o1n(int n, double k, Vec r, Vec theta, Vec phi){

  iVec *ThetaComp = new iVec(r.size());
  iVec *PhiComp   = new iVec(r.size());
  iVec *RComp     = new iVec(r.size());

  for (long unsigned int i = 0; i < r.size(); i++){

    double p = k*r[i];

    (*ThetaComp)[i] = (sin(phi[i]) * Taun(n, theta[i]) * ( Hankel(n, p) + p * Hankel_p(n, p)/p )  );

    (*PhiComp)[i]   = (  cos(phi[i]) * Pin(n, theta[i])  * ( Hankel(n, p) + p * Hankel_p(n, p)/p )  );

    (*RComp)[i]     = (    sin(phi[i]) * (complex128)( n * (n+1) ) * sin(theta[i]) * Pin(n, theta[i]) * Hankel(n, p)/p  );

  }

  return std::make_tuple(RComp, PhiComp, ThetaComp);
}


class VectorSphericalHarmonics3{

  public:

    int n;
    double k;


    VectorSphericalHarmonics3(int n, double Wavelength){
      this->n = n;
      this->k = 2 * PI /Wavelength;

    }


    complex128 ZFunc(std::size_t n, double x){ return Hankel(this->n, x);}

    complex128 ZFunc_p(std::size_t n, double x){ return Hankel_p(this->n, x);}


    std::tuple<iVec*, iVec*, iVec*>
    M_o1n(Vec r, Vec theta, Vec phi){
      iVec *ThetaComp = new iVec(r.size());
      iVec *PhiComp   = new iVec(r.size());
      iVec *RComp     = new iVec(r.size());

      for (long unsigned int i = 0; i < r.size(); i++){
        double p = this->k*r[i];
        ThetaComp->push_back(cos(phi[i]) * Pin(this->n, theta[i]) * this->ZFunc(this->n,p)  );
        PhiComp->push_back( -sin(phi[i]) * Taun(this->n, theta[i]) * this->ZFunc(this->n, p) );
        RComp->push_back(0.);
      }

      return std::make_tuple(RComp, PhiComp, ThetaComp);
    };

    std::tuple<iVec*, iVec*, iVec*>
    M_e1n(Vec r, Vec theta, Vec phi){
      iVec *ThetaComp = new iVec(r.size());
      iVec *PhiComp   = new iVec(r.size());
      iVec *RComp     = new iVec(r.size());

      for (long unsigned int i = 0; i < r.size(); i++){
        double p = this->k*r[i];
        ThetaComp->push_back(-sin(phi[i]) * Pin(this->n, theta[i]) * this->ZFunc(this->n, p) );
        PhiComp->push_back(  -cos(phi[i]) * Taun(this->n, theta[i]) * this->ZFunc(this->n, p) );
        RComp->push_back(0.);
      }

      return std::make_tuple(RComp, PhiComp, ThetaComp);
    };

    std::tuple<iVec*, iVec*, iVec*>
    N_o1n(Vec r, Vec theta, Vec phi){
      iVec *ThetaComp = new iVec(r.size());
      iVec *PhiComp   = new iVec(r.size());
      iVec *RComp     = new iVec(r.size());

      for (long unsigned int i = 0; i < r.size(); i++){
        double p = this->k*r[i];
        ThetaComp->push_back(sin(phi[i]) * Taun(this->n, theta[i]) * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p)/p )  );
        PhiComp->push_back(  cos(phi[i]) * Pin(this->n, theta[i])  * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p)/p )  );
        RComp->push_back(    sin(phi[i]) * this->n * (this->n+1) * sin(theta[i]) * Pin(this->n, theta[i]) * ZFunc(this->n, p)/p  );
      }

      return std::make_tuple(RComp, PhiComp, ThetaComp);
    };

    std::tuple<iVec*, iVec*, iVec*>
    N_e1n(Vec r, Vec theta, Vec phi){

      iVec *ThetaComp = new iVec(r.size());
      iVec *PhiComp   = new iVec(r.size());
      iVec *RComp     = new iVec(r.size());

      for (long unsigned int i = 0; i < r.size(); i++){

        double p = this->k*r[i];


        (*ThetaComp)[i] =  cos(phi[i]) * Taun(this->n, theta[i]) * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p) ) / p ;

        (*PhiComp)[i]   = -sin(phi[i]) * Pin(this->n, theta[i])  * ( this->ZFunc(this->n, p) + p * this->ZFunc_p(this->n, p) ) / p ;

        (*RComp)[i]     =  cos(phi[i]) * (complex128)this->n * (complex128)(this->n+1) * sin(theta[i]) * Pin(this->n, theta[i]) * ZFunc(this->n, p) / p ;

      }

      return std::make_tuple(RComp, PhiComp, ThetaComp);
    };

};

// -
