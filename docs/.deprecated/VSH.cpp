

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
