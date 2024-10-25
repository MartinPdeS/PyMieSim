#pragma once

#include "single/headers/fibonacci_mesh.h"

#define PI (double)3.14159265358979323846264338

class FullSteradian : public BaseMesh
{
    public:
        double dTheta, dPhi;

        FullSteradian(const size_t sampling, const double radius = 1.0) : BaseMesh(sampling, radius)
        {
            dTheta = 2.0 * PI / (sampling-1);
            dPhi   = 1.0 * PI / (sampling-1);

            for (size_t p=0; p<sampling; p++)
                spherical_coordinates.phi.push_back(p * dPhi - PI / 2.0);

            for (size_t t=0; t<sampling; t++)
                spherical_coordinates.theta.push_back(t * dTheta - PI / 1.0);
        }

        template<typename dtype>
        dtype get_integral(std::vector<dtype>& vector) const
        {
            dtype integral = 0;

            for (size_t p=0; p<sampling; p++)
                for (size_t t=0; t<sampling; t++)
                    integral += vector[p*sampling + t] * sin(spherical_coordinates.phi[p] + PI/2.0) * dPhi * dTheta;

            return integral;
        }


        template<typename dtype>
        dtype get_cos_integral(const std::vector<dtype>& vector) const
        {
            dtype integral = 0;

            for (size_t p=0; p<sampling; p++)
                for (size_t t=0; t<sampling; t++)
                    integral += vector[p*sampling + t] * cos(spherical_coordinates.phi[p] + PI/2.0) * sin(spherical_coordinates.phi[p] + PI/2.0) * dPhi * dTheta;

            return integral;
        }

        double get_integral() const
        {
            double integral = 0;

            for (auto phi : spherical_coordinates.phi)
                for (size_t theta_idx = 0; theta_idx < spherical_coordinates.theta.size(); ++theta_idx)
                    integral += sin(phi+PI/2.0) * dPhi * dTheta;

            return integral;
        }

};

// --
