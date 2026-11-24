#include "full_mesh.h"


FullSteradian::FullSteradian(const size_t sampling, const double radius)
 : sampling(sampling), radius(radius)
{
    total_size = sampling * sampling;
    dTheta = 2.0 * Constants::PI / (sampling - 1);
    dPhi   = 1.0 * Constants::PI / (sampling - 1);

    for (size_t p=0; p<sampling; p++)
        this->spherical.phi.push_back(p * dPhi - Constants::PI / 2.0);
    for (size_t t=0; t<sampling; t++)
        this->spherical.theta.push_back(t * dTheta - Constants::PI / 1.0);

    this->spherical.shape = {sampling,};

    for (size_t p=0; p<sampling; p++)
        for (size_t t=0; t<sampling; t++)
        {
            double phi = this->spherical.phi[p];
            double theta = this->spherical.theta[t];

            this->spherical_mesh.r.push_back(radius);
            this->spherical_mesh.phi.push_back(phi);
            this->spherical_mesh.theta.push_back(theta);
            this->spherical_mesh.shape = {sampling, sampling};

            this->cartesian_mesh.x.push_back(radius * sin(phi) * cos(theta));
            this->cartesian_mesh.y.push_back(radius * sin(phi) * sin(theta));
            this->cartesian_mesh.z.push_back(radius * cos(phi));
            this->cartesian_mesh.shape = {sampling, sampling};
        }
}

template <typename dtype>
dtype FullSteradian::get_integral(std::vector<dtype> &vector) const {
    dtype integral = 0;

    for (size_t elevation_index = 0; elevation_index < sampling; ++elevation_index) {
        const double weight = std::cos(spherical.phi[elevation_index]) * dPhi * dTheta;
        // The endpoints of the azimuth grid represent the same direction.
        for (size_t azimuth_index = 0; azimuth_index + 1 < sampling; ++azimuth_index) {
            integral += vector[elevation_index * sampling + azimuth_index] * weight;
        }
    }

    return integral;
}

template <typename dtype>
dtype FullSteradian::get_cos_integral(const std::vector<dtype> &vector) const {
    dtype integral = 0;

    for (size_t elevation_index = 0; elevation_index < sampling; ++elevation_index) {
        const double elevation = spherical.phi[elevation_index];
        // The cosine of the scattering angle from +z equals the sine of elevation.
        const double weight = std::sin(elevation) * std::cos(elevation) * dPhi * dTheta;
        for (size_t azimuth_index = 0; azimuth_index + 1 < sampling; ++azimuth_index) {
            integral += vector[elevation_index * sampling + azimuth_index] * weight;
        }
    }

    return integral;
}

double FullSteradian::get_integral() const {
    double integral = 0.0;

    for (const double elevation : spherical.phi) {
        integral += std::cos(elevation) * dPhi;
    }

    return 2.0 * Constants::PI * integral;
}

// Explicit instantiations for the types you want:
template double FullSteradian::get_integral<double>(std::vector<double>&) const;
template complex128 FullSteradian::get_integral<complex128>(std::vector<complex128>&) const;

template double FullSteradian::get_cos_integral<double>(const std::vector<double>&) const;
template complex128 FullSteradian::get_cos_integral<complex128>(const std::vector<complex128>&) const;
// --
