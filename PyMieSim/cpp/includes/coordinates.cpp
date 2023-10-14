#ifndef COORDINATE_H
#define COORDINATE_H

#include "definitions.cpp"
#include "special_function.cpp"
#include "numpy_interface.cpp"


  struct SCoordinate { double R, Phi, Theta; };


  struct VectorField
  {
    size_t Sampling, Size;
    IVector Shape;
    DVector Data;

    VectorField(){}
    VectorField(const DVector& Vector) : Sampling( 1 ), Size( 3 ), Shape( { 1, 3 } ), Data( Vector ) {}

    VectorField(const size_t &Sampling) : Sampling( Sampling ), Size( 3*Sampling ), Shape( { Sampling, 3 } ), Data( DVector( 3*Sampling ) ) {}
    ndarray GetNumpy(){ return vector_to_ndarray_copy((*this).Data, (*this).Shape); }

    double      &operator[](size_t i) { return Data[i]; }
    double      &operator()(size_t& i, size_t j) { return Data[i*3 + j]; }
    double      &operator()(size_t&& i, size_t&& j) { return Data[i*3 + j]; }

    VectorField operator+(VectorField& Other)
    {
      VectorField Output(Sampling);

      for (size_t i=0; i<3*Sampling; i++)
          Output[i] = (*this)[i] + Other[i];

      return Output;
    }

    DVector ScalarProduct(DVector& Base)
    {
      DVector Output(Sampling);

      double x0  = Base[0],
             y0  = Base[1],
             z0  = Base[2];


      for (size_t i=0; i<Sampling; i++)
      {
        double x  = (*this)(i,0),
               y  = (*this)(i,1),
               z  = (*this)(i,2);

        double Proj = x * x0 + y * y0 + z * z0;

        Output[i] = Proj;
      }

      return Output;
    }


    DVector ScalarProduct(VectorField& Base)
    {
      DVector Output(Sampling);

      double x0  = Base(0,0),
             y0  = Base(0,1),
             z0  = Base(0,2);


      for (size_t i=0; i<Sampling; i++)
      {
        double x  = (*this)(i,0),
               y  = (*this)(i,1),
               z  = (*this)(i,2);

        double Proj = x * x0 + y * y0 + z * z0;

        Output[i] = Proj;
      }

      return Output;
    }


    VectorField Projection(DVector& Base)
    {
      VectorField Output(Sampling);

      double x0 = Base[0],
             y0 = Base[1],
             z0 = Base[2];


      for (size_t i=0; i<Sampling; i++)
      {
        double x  = (*this)(i,0),
               y  = (*this)(i,1),
               z  = (*this)(i,2);

        double Proj = x * x0 + y * y0 + z * z0;
        double Norm = x0 * x0 + y0 * y0 + z0 * z0;

        Output(i,0) = Proj/Norm * x0;
        Output(i,1) = Proj/Norm * y0;
        Output(i,2) = Proj/Norm * z0;
      }
      return Output;
    }

    VectorField Projection(VectorField& Base)
    {
      VectorField Output(Sampling);

      double x0 = Base(0,0),
             y0 = Base(0,1),
             z0 = Base(0,2);
      double Norm = x0 * x0 + y0 * y0 + z0 * z0;

      for (size_t i=0; i<Sampling; i++)
      {
        double x  = (*this)(i,0),
               y  = (*this)(i,1),
               z  = (*this)(i,2);

        double Proj = x * x0 + y * y0 + z * z0;

        Output(i,0) = Proj/Norm * x0;
        Output(i,1) = Proj/Norm * y0;
        Output(i,2) = Proj/Norm * z0;
      }
      return Output;
    }


    VectorField Projection(VectorField& Base0, VectorField& Base1)
    {
      VectorField Output(Sampling);

      double x0 = Base0(0,0),
             y0 = Base0(0,1),
             z0 = Base0(0,2),
             x1 = Base1(0,0),
             y1 = Base1(0,1),
             z1 = Base1(0,2);

      double Norm0 = x0 * x0 + y0 * y0 + z0 * z0;
      double Norm1 = x1 * x1 + y1 * y1 + z1 * z1;

      for (size_t i=0; i<Sampling; i++)
      {
        double xp  = (*this)(i,0),
               yp  = (*this)(i,1),
               zp  = (*this)(i,2);

        double Proj0 = xp * x0 + yp * y0 + zp * z0,
               Proj1 = xp * x1 + yp * y1 + zp * z1;

        Output(i,0) = Proj0/Norm0 * x0 + Proj1/Norm1 * x1;
        Output(i,1) = Proj0/Norm0 * y0 + Proj1/Norm1 * y1;
        Output(i,2) = Proj0/Norm0 * z0 + Proj1/Norm1 * z1;
      }
      return Output;
    }


    void Normalize()
    {
      for (size_t i=0; i<Sampling; i++)
      {
        double x = (*this)(i,0), y = (*this)(i,1), z = (*this)(i,2);
        double norm = sqrt( pow(x, 2) + pow(y, 2) + pow(z, 2) );
        if (norm != 0.0)
        {
          (*this)(i,0) /= norm;
          (*this)(i,1) /= norm;
          (*this)(i,2) /= norm;
        }
      }
    }

    void mx_rot_y(double theta)
    {
      Matrix3 M
      {
        { cos(theta),  0.,    sin(theta) },
        { 0.,          1.,    0.         },
        {-sin(theta),  0.,    cos(theta) }
      };
    mx_apply(M);
    }

    void mx_rot_x(double gamma)
    {
      Matrix3 M
      {
        {1.,          0.,            0.         },
        {0.,          cos(gamma),   -sin(gamma) },
        {0.,          sin(gamma),    cos(gamma) }
      };
    mx_apply(M);
    }

    void mx_rot_z(double phi)
    {
      Matrix3 M
      {
        {cos(phi),    -sin(phi),   0.     },
        {sin(phi),     cos(phi),   0.     },
        {0.,             0.,       1.     }
      };
    mx_apply(M);
    }

    void mx_apply(Matrix3 M)
    {
      double tempx, tempy, tempz;
      for (size_t i = 0; i < Sampling; i++){

          tempx = M[0][0] * (*this)(i,0) + M[0][1] * (*this)(i,1) + M[0][2] * (*this)(i,2);
          tempy = M[1][0] * (*this)(i,0) + M[1][1] * (*this)(i,1) + M[1][2] * (*this)(i,2);
          tempz = M[2][0] * (*this)(i,0) + M[2][1] * (*this)(i,1) + M[2][2] * (*this)(i,2);

          (*this)(i,0) = tempx;
          (*this)(i,1) = tempy;
          (*this)(i,2) = tempz;
        }
    }

  };




  struct Spherical
  {
    DVector Phi, Theta, R;
    size_t Sampling;

    Spherical(){}
    Spherical(const size_t &Sampling)
                     : Phi(DVector(Sampling)),
                       Theta(DVector(Sampling)),
                       R(DVector(Sampling)),
                       Sampling(Sampling)
                       {}
    ndarray get_r_py()    {return vector_to_ndarray_copy(R);}
    ndarray get_phi_py()  {return vector_to_ndarray_copy(Phi);}
    ndarray get_theta_py(){return vector_to_ndarray_copy(Theta);}
  };


  struct Cartesian
  {
    DVector X, Y, Z;
    size_t Sampling;

    Cartesian(){}
    Cartesian(const size_t &Sampling)
                     : X(DVector(Sampling)),
                       Y(DVector(Sampling)),
                       Z(DVector(Sampling)),
                       Sampling(Sampling)
                       {}

    ndarray get_x_py(){return vector_to_ndarray_copy(X);}
    ndarray get_y_py(){return vector_to_ndarray_copy(Y);}
    ndarray get_z_py(){return vector_to_ndarray_copy(Z);}

    Spherical Cart2Sph()
    {
      Spherical Output(Sampling);

      for (size_t i = 0; i < Sampling; i++){
         Output.R[i] = sqrt( pow( X[i],2) + pow(Y[i],2) + pow(Z[i],2) ) ;
         Output.Theta[i] = atan2( Y[i],  X[i] );
         Output.Phi[i] = asin( Z[i] / Output.R[i] );
         }
        return Output;
    }

    void mx_rot_y(double theta)
    {
      Matrix3 M
      {
        { cos(theta),  0.,    sin(theta) },
        { 0.,          1.,    0.         },
        {-sin(theta),  0.,    cos(theta) }
      };
    mx_apply(M);
    }

    void mx_rot_x(double gamma)
    {
      Matrix3 M
      {
        {1.,          0.,            0.         },
        {0.,          cos(gamma),   -sin(gamma) },
        {0.,          sin(gamma),    cos(gamma) }
      };
    mx_apply(M);
    }

    void mx_rot_z(double phi)
    {
      Matrix3 M
      {
        {cos(phi),    -sin(phi),   0.     },
        {sin(phi),     cos(phi),   0.     },
        {0.,             0.,       1.     }
      };
    mx_apply(M);
    }

    void mx_apply(Matrix3 M)
    {
      double tempx, tempy, tempz;
      for (size_t i = 0; i < Sampling; i++){

          tempx = M[0][0] * X[i] + M[0][1] * Y[i] + M[0][2] * Z[i];
          tempy = M[1][0] * X[i] + M[1][1] * Y[i] + M[1][2] * Z[i];
          tempz = M[2][0] * X[i] + M[2][1] * Y[i] + M[2][2] * Z[i];

          X[i] = tempx;
          Y[i] = tempy;
          Z[i] = tempz;
        }
    }

  };









#endif


// --
