#ifndef COORDINATE_H
#define COORDINATE_H

#include "definitions.cpp"
#include "special_function.cpp"
#include "numpy_interface.cpp"


  struct SCoordinate { double R, Phi, Theta; };


  struct VectorField
  {
    size_t sampling, Size;
    IVector shape;
    DVector Data;

    VectorField(){}
    VectorField(const DVector& Vector) : sampling(1), Size(3), shape({1, 3}), Data(Vector) {}

    VectorField(const size_t &sampling) : sampling(sampling), Size(3*sampling), shape( {sampling, 3} ), Data( DVector(3 * sampling) ) {}
    ndarray get_numpy(){ return vector_to_ndarray_copy((*this).Data, (*this).shape); }

    double &operator[](size_t i) { return Data[i]; }
    double &operator()(size_t& i, size_t j) { return Data[i*3 + j]; }
    double &operator()(size_t&& i, size_t&& j) { return Data[i*3 + j]; }

    VectorField operator+(VectorField& Other)
    {
      VectorField output(sampling);

      for (size_t i=0; i<3*sampling; i++)
          output[i] = (*this)[i] + Other[i];

      return output;
    }

    DVector get_scalar_product(DVector& Base)
    {
      DVector output(sampling);

      double
        x0  = Base[0],
        y0  = Base[1],
        z0  = Base[2];


      for (size_t i=0; i<sampling; i++)
      {
        double
          x  = (*this)(i,0),
          y  = (*this)(i,1),
          z  = (*this)(i,2),
          projection = x * x0 + y * y0 + z * z0;

        output[i] = projection;
      }

      return output;
    }


    std::vector<double> get_scalar_product(VectorField& base_vector)
    {
      std::vector<double>
        scalar_product(sampling);

      double
        x0 = base_vector(0,0),
        y0 = base_vector(0,1),
        z0 = base_vector(0,2);


      for (size_t i=0; i<sampling; i++)
      {
        double
          x = (*this)(i,0),
          y = (*this)(i,1),
          z = (*this)(i,2),
          projection = x * x0 + y * y0 + z * z0;

        scalar_product[i] = projection;
      }

      return scalar_product;
    }


    VectorField get_projection(const std::vector<double>& base_vector)
    {
      VectorField output(sampling);

      double
        x0 = base_vector[0],
        y0 = base_vector[1],
        z0 = base_vector[2];


      for (size_t i=0; i<sampling; i++)
      {
        double
          x = (*this)(i,0),
          y = (*this)(i,1),
          z = (*this)(i,2),
          projection = x * x0 + y * y0 + z * z0,
          norm = x0 * x0 + y0 * y0 + z0 * z0;

        output(i,0) = projection / norm * x0;
        output(i,1) = projection / norm * y0;
        output(i,2) = projection / norm * z0;
      }
      return output;
    }

    VectorField get_projection(VectorField& base_vector)
    {
      VectorField output(sampling);

      double
        x0 = base_vector(0,0),
        y0 = base_vector(0,1),
        z0 = base_vector(0,2),
        norm = x0 * x0 + y0 * y0 + z0 * z0;

      for (size_t i=0; i<sampling; i++)
      {
        double
          x  = (*this)(i,0),
          y  = (*this)(i,1),
          z  = (*this)(i,2),
          projection = x * x0 + y * y0 + z * z0;

        output(i,0) = projection / norm * x0;
        output(i,1) = projection / norm * y0;
        output(i,2) = projection / norm * z0;
      }
      return output;
    }


    VectorField get_projection(VectorField& base_vector_0, VectorField& base_vector_1)
    {
      VectorField output(sampling);

      double
        x0 = base_vector_0(0,0),
        y0 = base_vector_0(0,1),
        z0 = base_vector_0(0,2),
        x1 = base_vector_1(0,0),
        y1 = base_vector_1(0,1),
        z1 = base_vector_1(0,2),
        norm_0 = x0 * x0 + y0 * y0 + z0 * z0,
        norm_1 = x1 * x1 + y1 * y1 + z1 * z1;

      for (size_t i=0; i<sampling; i++)
      {
        double
          xp  = (*this)(i,0),
          yp  = (*this)(i,1),
          zp  = (*this)(i,2),
          projection_0 = xp * x0 + yp * y0 + zp * z0,
          projection_1 = xp * x1 + yp * y1 + zp * z1,
          factor_0 = projection_0 / norm_0,
          factor_1 = projection_1 / norm_1;

        output(i,0) = factor_0 * x0 + factor_1 * x1;
        output(i,1) = factor_0 * y0 + factor_1 * y1;
        output(i,2) = factor_0 * z0 + factor_1 * z1;
      }
      return output;
    }


    void normalize()
    {
      for (size_t i=0; i<sampling; i++)
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
      double
        tempx,
        tempy,
        tempz;

      for (size_t i = 0; i < sampling; i++){

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
    std::vector<double>
      Phi,
      Theta,
      R;

    size_t sampling;

    Spherical(){}
    Spherical(const size_t &sampling)
    : Phi(std::vector<double>(sampling)),
      Theta(std::vector<double>(sampling)),
      R(std::vector<double>(sampling)),
      sampling(sampling)
      {}

    ndarray get_r_py() {return vector_to_ndarray_copy(R);}
    ndarray get_phi_py(){return vector_to_ndarray_copy(Phi);}
    ndarray get_theta_py(){return vector_to_ndarray_copy(Theta);}
  };


  struct Cartesian
  {
    std::vector<double>
      X,
      Y,
      Z;

    size_t sampling;

    Cartesian(){}
    Cartesian(const size_t &sampling)
    : X(std::vector<double>(sampling)),
      Y(std::vector<double>(sampling)),
      Z(std::vector<double>(sampling)),
      sampling(sampling)
      {}

    void set_x_py(const std::vector<double> &value){X = value;}
    void set_y_py(const std::vector<double> &value){Y = value;}
    void set_z_py(const std::vector<double> &value){Z = value;}

    ndarray get_x_py(){return vector_to_ndarray_copy(X);}
    ndarray get_y_py(){return vector_to_ndarray_copy(Y);}
    ndarray get_z_py(){return vector_to_ndarray_copy(Z);}

    Spherical cartesian_to_spherical()
    {
      Spherical output(sampling);

      for (size_t i = 0; i < sampling; i++){
         output.R[i] = sqrt( pow( X[i],2) + pow(Y[i],2) + pow(Z[i],2) ) ;
         output.Theta[i] = atan2( Y[i],  X[i] );
         output.Phi[i] = asin( Z[i] / output.R[i] );
         }
        return output;
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
      for (size_t i = 0; i < sampling; i++){

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
