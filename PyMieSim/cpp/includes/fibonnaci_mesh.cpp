#ifndef FIBONACCI_H
#define FIBONACCI_H

  #include "coordinates.cpp"
  #include "numpy_interface.cpp"


  class FibonacciMesh
  {

    public:
        size_t Sampling, TrueSample;
        Spherical SCoord;
        Cartesian CCoord, BaseCCoord;
        VectorField Perp, Para;

        std::vector<double>
            HPara,
            VPara,
            HPerp,
            VPerp;

        double
            MaxAngle,
            phi_offset,
            gamma_offset,
            dOmega,
            Omega;

        VectorField
            VVector = VectorField({1,0,0}),
            HVector = VectorField({0,1,0});

        void
            mx_rot_y(double theta),
            mx_rot_x(double theta),
            mx_rot_z(double theta),
            mx_apply(Matrix3 M),
            GenCoordinate(),
            Cart2Sph();

        ndarray GetPara() {return Para.GetNumpy();};
        ndarray GetPerp() {return Perp.GetNumpy();};

        ndarray get_H_Para() {return vector_to_ndarray_copy(HPara);};
        ndarray get_V_Para() {return vector_to_ndarray_copy(VPara);};
        ndarray get_H_Perp() {return vector_to_ndarray_copy(HPerp);};
        ndarray get_V_Perp() {return vector_to_ndarray_copy(VPerp);};

        ndarray get_x_py() {return CCoord.get_x_py();};
        ndarray get_y_py() {return CCoord.get_y_py();};
        ndarray get_z_py() {return CCoord.get_z_py();};


        ndarray get_r_py() {return SCoord.get_r_py();};
        ndarray get_phi_py() {return SCoord.get_phi_py();};
        ndarray get_theta_py() { return SCoord.get_theta_py();};

        FibonacciMesh(){}

        FibonacciMesh(const int    &Sampling,
                      const double &MaxAngle,
                      const double &phi_offset,
                      const double &gamma_offset)
                      : Sampling(Sampling),
                        MaxAngle(MaxAngle),
                        phi_offset(phi_offset),
                        gamma_offset(gamma_offset)
        {
        CCoord = Cartesian(Sampling);
        compute_mesh();
        BaseCCoord = CCoord;
        Rotations();
        SCoord = CCoord.Cart2Sph() ;
        compute_vector_field();
        compute_HV_projection();
        }


  void Rotations()
  {
    if (gamma_offset != 0.)
    {
        CCoord.mx_rot_x(gamma_offset);
        VVector.mx_rot_x(gamma_offset);
        HVector.mx_rot_x(gamma_offset);
    }

    if (phi_offset   != 0.)
    {
        CCoord.mx_rot_y(phi_offset);
        VVector.mx_rot_y(phi_offset);
        HVector.mx_rot_y(phi_offset);
    }
  }


  void compute_vector_field()
  {
    Para = VectorField(Sampling);
    Perp = VectorField(Sampling);

    for (size_t i = 0; i < Sampling; i++)
    {
        Perp(i,0) = -CCoord.Y[i] ;
        Perp(i,1) = CCoord.X[i] ;
        Perp(i,2) = 0.0;

        Para(i,0) = CCoord.X[i] * CCoord.Z[i] ;
        Para(i,1) = CCoord.Y[i] * CCoord.Z[i] ;
        Para(i,2) = -( pow( CCoord.X[i], 2 ) + pow( CCoord.Y[i], 2 ) );
    }
    Para.Normalize();
    Perp.Normalize();
  }




    void compute_HV_projection()
    {
        HPara = Para.ScalarProduct(HVector);
        VPara = Para.ScalarProduct(VVector);

        HPerp = Perp.ScalarProduct(HVector);
        VPerp = Perp.ScalarProduct(VVector);
    }


    double NA2Angle(const double &NA)
    {
        if (NA <= 1.0)
            return asin(NA);

        if (NA >= 1.0)
            return asin(NA-1.0) + PI/2.0;

        return 1.0;
    }


  void compute_mesh()
  {
    this->get_properties();

    double phi  = PI * (3. - sqrt(5.));  // golden angle = 2.39996322972865332

    for (size_t i = 0; i < TrueSample; i++){
        double theta    = phi * i;
        CCoord.Z[i]     = 1 - (i / (double)( TrueSample - 1) ) * 2 ;
        double radius   = sqrt(1 - CCoord.Z[i] * CCoord.Z[i]);
        CCoord.X[i]     = cos(theta) * radius ;
        CCoord.Y[i]     = sin(theta) * radius ;

        if (i == (size_t)Sampling-1)
            break;
    }
  }


  void get_properties()
  {
    double
        solid_angle = abs( 2. * PI * ( cos(MaxAngle) - 1. ) ),   //cos(0) =1
        ratio = ( 4. * PI / solid_angle );

    size_t _TrueSample = (size_t) ( Sampling * ratio );

    dOmega = 4. * PI / _TrueSample;
    Omega = dOmega * Sampling ;
    TrueSample = _TrueSample;
  }

};



class FullSteradian
{
    public:
    size_t Sampling;

    Spherical SCoord;

    double dTheta, dPhi;

    FullSteradian(size_t Sampling) : Sampling(Sampling), SCoord(Sampling)
    {
        dTheta = 2.0 * PI / (Sampling-1);
        dPhi   = 1.0 * PI / (Sampling-1);

        for (size_t p=0; p<Sampling; p++)
            SCoord.Phi[p]   = p * dPhi   - PI/2.0;

        for (size_t t=0; t<Sampling; t++)
            SCoord.Theta[t] = t * dTheta - PI/1.0;
    }

    template<typename T> T Integral(std::vector<T>& Vector)
    {
        T integral = 0;

        for (size_t p=0; p<Sampling; p++)
            for (size_t t=0; t<Sampling; t++)
                integral += Vector[p*Sampling + t] * sin(SCoord.Phi[p] + PI/2.0) * dPhi * dTheta;

        return integral;
    }


    template<typename T> T IntegralCos(std::vector<T>& Vector)
    {
        T integral = 0;

        for (size_t p=0; p<Sampling; p++)
            for (size_t t=0; t<Sampling; t++)
                integral += Vector[p*Sampling + t] * cos(SCoord.Phi[p] + PI/2.0) * sin(SCoord.Phi[p] + PI/2.0) * dPhi * dTheta;

        return integral;
    }

    double Integral()
    {
        double integral = 0;

        for (auto phi : SCoord.Phi)
            for (auto theta : SCoord.Theta)
                integral += sin(phi+PI/2.0) * dPhi * dTheta;

        return integral;
    }
};


#endif

// --
