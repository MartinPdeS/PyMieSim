#include <iostream>
#include <vector>
#include <complex>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

typedef std::complex<double> complex128;
typedef std::vector<double> Vec;
typedef py::array_t<double> ndarray;
typedef std::vector<std::vector<double>> Matrix3;

#define J complex128(0.0,1.0)
#define PI (double) 3.14159265358979323846264338

using namespace std;



class FibonacciMesh{

  public:
    int      Samples, TrueSample;

    ndarray  x, y, z, R, Phi, Theta, Sin, PhiBase, ThetaBase;

    double   MaxAngle,
             PhiOffset,
             GammaOffset,
            *x_data,
            *y_data,
            *z_data,
            *R_data,
            *Phi_data,
            *Theta_data,
            *PhiBase_data,
            *ThetaBase_data,
            *Sin_data,
             dOmega,
             Omega;

    void     GenMesh(),
             mx_rot_y(double theta),
             mx_rot_x(double theta),
             mx_rot_z(double theta),
             mx_apply(Matrix3 M),
             GetBaseCoordinate(),
             GetProperties(),
             GenCoordinate(),
             Cart2Sph();

    ndarray& Getx()           { return this->x;         };
    ndarray& Gety()           { return this->y;         };
    ndarray& Getz()           { return this->z;         };
    ndarray& GetR()           { return this->R;         };
    ndarray& GetPhi()         { return this->Phi;       };
    ndarray& GetTheta()       { return this->Theta;     };
    ndarray& GetPhiBase()     { return this->PhiBase;   };
    ndarray& GetThetaBase()   { return this->ThetaBase; };

    ndarray& GetSin()         { return this->Sin;       };

    double&  GetdOmega()      { return this->dOmega;    };
    double&  GetOmega()       { return this->Omega;     };

    FibonacciMesh(int    Samples,
                  double MaxAngle,
                  double PhiOffset,
                  double GammaOffset)
    {
      this->Samples     = Samples;
      this->MaxAngle    = MaxAngle;
      this->PhiOffset   = PhiOffset;
      this->GammaOffset = GammaOffset;

      this->GenCoordinate();

      this->GenMesh();

      this->Cart2Sph();

      this->GetBaseCoordinate();

      if (this->GammaOffset != 0.) {this->mx_rot_x(this->GammaOffset); }

      if (this->PhiOffset   != 0.) {this->mx_rot_y(this->PhiOffset);   }

      this->Cart2Sph();
    }
};


void
FibonacciMesh::Cart2Sph(){

  for (auto i = 0; i < this->Samples; i++){
      double rho = sqrt( pow(this->z_data[i],2) + pow(this->x_data[i],2) ) ;

      this->R_data[i]       = sqrt( pow( this->x_data[i],2) + pow(this->y_data[i],2) + pow(this->z_data[i],2) ) ;
      this->Theta_data[i]   = atan2( this->y_data[i],  this->x_data[i] );
      this->Phi_data[i]     = asin( this->z_data[i] / this->R_data[i] );
      this->Sin_data[i]     = abs( sin(Phi_data[i]) );
      }
}


void
FibonacciMesh::GenMesh()
{
    this->GetProperties();

    double phi  = PI * (3. - sqrt(5.));  // golden angle = 2.39996322972865332

    for (auto i = 0; i < this->TrueSample; i++){

       double theta  = phi * i;

       this->z_data[i] = 1 - (i / (double)( this->TrueSample - 1) ) * 2 ;

       double radius   = sqrt(1 - this->z_data[i] * this->z_data[i]);

       this->x_data[i] = cos(theta) * radius ;
       this->y_data[i] = sin(theta) * radius ;

       if (i == (int)this->Samples-1) break;
      }
}


void
FibonacciMesh::GetBaseCoordinate(){
  double temp;

  for (auto i = 0; i < this->Samples; i++){
      this->PhiBase_data[i]   = this->Phi_data[i];
      this->ThetaBase_data[i] = this->Theta_data[i];
      }
}


void
FibonacciMesh::GetProperties(){
  double SolidAngle   = abs( 2. * PI * ( cos(this->MaxAngle) - 1. ) ),   //cos(0) =1
         ratio        = ( 4. * PI / SolidAngle ),
         TrueSample   = (int) ( this->Samples * ratio );

  this->dOmega     = 4. * PI / TrueSample;
  this->Omega      = this->dOmega * this->Samples ;
  this->TrueSample = TrueSample;
}


void
FibonacciMesh::mx_rot_y(double theta){
  Matrix3 M
  {
    { cos(theta),  0.,    sin(theta) },
    { 0.,          1.,    0.         },
    {-sin(theta),  0.,    cos(theta) }
  };
this->mx_apply(M);
}


void
FibonacciMesh::mx_rot_x(double gamma){
  Matrix3 M
  {
    {1.,          0.,            0.         },
    {0.,          cos(gamma),   -sin(gamma) },
    {0.,          sin(gamma),    cos(gamma) }
  };
this->mx_apply(M);
}


void
FibonacciMesh::mx_rot_z(double phi){
  Matrix3 M
  {
    {cos(phi),    -sin(phi),   0.     },
    {sin(phi),     cos(phi),   0.     },
    {0.,             0.,       1.     }
  };
this->mx_apply(M);
}


void
FibonacciMesh::mx_apply(Matrix3 M){
  double tempx, tempy, tempz;
  for (auto i = 0; i < this->Samples; i++){

      tempx = M[0][0] * this->x_data[i] + M[0][1] * this->y_data[i] + M[0][2] * this->z_data[i];
      tempy = M[1][0] * this->x_data[i] + M[1][1] * this->y_data[i] + M[1][2] * this->z_data[i];
      tempz = M[2][0] * this->x_data[i] + M[2][1] * this->y_data[i] + M[2][2] * this->z_data[i];

      this->x_data[i] = tempx;
      this->y_data[i] = tempy;
      this->z_data[i] = tempz;
    }
}



void
FibonacciMesh::GenCoordinate()
{

  this->x              = ndarray(Samples);
  this->y              = ndarray(Samples);
  this->z              = ndarray(Samples);

  this->R              = ndarray(Samples);
  this->Phi            = ndarray(Samples);
  this->Theta          = ndarray(Samples);
  this->Sin            = ndarray(Samples);

  this->PhiBase        = ndarray(Samples);
  this->ThetaBase      = ndarray(Samples);

  this->x_data         = this->x.mutable_data();
  this->y_data         = this->y.mutable_data();
  this->z_data         = this->z.mutable_data();

  this->R_data         = this->R.mutable_data();
  this->Phi_data       = this->Phi.mutable_data();
  this->Theta_data     = this->Theta.mutable_data();
  this->Sin_data       = this->Sin.mutable_data();

  this->PhiBase_data   = this->PhiBase.mutable_data();
  this->ThetaBase_data = this->ThetaBase.mutable_data();
}


PYBIND11_MODULE(_utils, module) {
    module.doc() = "LGeneralized Lorenz-Mie Theory (GLMT) c++ binding module for light scattering from a spherical scatterer";


      py::class_<FibonacciMesh>(module, "FibonacciMeshCpp")
      .def(py::init<int, double, double, double>())
      .def_property("x",          &FibonacciMesh::Getx,         &FibonacciMesh::Getx)
      .def_property("y",          &FibonacciMesh::Gety,         &FibonacciMesh::Gety)
      .def_property("z",          &FibonacciMesh::Getz,         &FibonacciMesh::Getz)
      .def_property("r",          &FibonacciMesh::GetR,         &FibonacciMesh::Getx)
      .def_property("phi",        &FibonacciMesh::GetPhi,       &FibonacciMesh::Gety)
      .def_property("theta",      &FibonacciMesh::GetTheta,     &FibonacciMesh::Getz)
      .def_property("PhiBase",    &FibonacciMesh::GetPhiBase,   &FibonacciMesh::Getz)
      .def_property("ThetaBase",  &FibonacciMesh::GetThetaBase, &FibonacciMesh::Getz)
      .def_property("SinMesh",    &FibonacciMesh::GetSin,       &FibonacciMesh::Getz)
      .def_property("dOmega",     &FibonacciMesh::GetdOmega,    &FibonacciMesh::Getz)
      .def_property("Omega",      &FibonacciMesh::GetOmega,     &FibonacciMesh::Getz);
}





















// --
