#pragma once

#include <vector>
#include <complex>
#include "utils/base_mesh.h"

#define PI (double)3.14159265358979323846264338

typedef std::complex<double> complex128;

class FullSteradian : public BaseMesh
{
    public:
        double dTheta, dPhi;
        Cartesian cartesian;
        Spherical spherical;

        /**
         * @brief Default constructor for FullSteradian mesh.
         * Initializes the mesh with a specified sampling size and radius.
         * @param sampling The number of samples to reserve for the FullSteradian mesh.
         * @param radius The radius of the mesh (default is 1.0).
         */
        FullSteradian(const size_t sampling, const double radius = 1.0);

        /**
         * @brief Computes the integral of a vector field over the full steradian.
         * @param vector The vector field to integrate.
         * @return The computed integral value.
         */
        template<typename dtype> dtype get_integral(std::vector<dtype>& vector) const;

        /**
         * @brief Computes the integral of a vector field over the full steradian.
         * @param vector The vector field to integrate.
         * @return The computed integral value.
         */
        template<typename dtype> dtype get_cos_integral(const std::vector<dtype>& vector) const;

        /**
         * @brief Computes the integral over the full steradian.
         * @return The computed integral value.
         */
        double get_integral() const;

};

// --
