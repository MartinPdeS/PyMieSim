#include "./full_mesh.h"


FullSteradian::FullSteradian(const size_t sampling, const double radius)
 : sampling(sampling), radius(radius)
{
    dTheta = 2.0 * PI / (sampling - 1);
    dPhi   = 1.0 * PI / (sampling - 1);

    for (size_t p=0; p<sampling; p++)
        this->spherical.phi.push_back(p * dPhi - PI / 2.0);

    for (size_t t=0; t<sampling; t++)
        this->spherical.theta.push_back(t * dTheta - PI / 1.0);
}

template<typename dtype>
dtype FullSteradian::get_integral(std::vector<dtype>& vector) const
{
    dtype integral = 0;

    for (size_t p=0; p<sampling; p++)
        for (size_t t=0; t<sampling; t++)
            integral += vector[p*sampling + t] * sin(this->spherical.phi[p] + PI/2.0) * dPhi * dTheta;

    return integral;
}


template<typename dtype>
dtype FullSteradian::get_cos_integral(const std::vector<dtype>& vector) const
{
    dtype integral = 0;

    for (size_t p=0; p<sampling; p++)
        for (size_t t=0; t<sampling; t++)
            integral += vector[p*sampling + t] * cos(this->spherical.phi[p] + PI/2.0) * sin(this->spherical.phi[p] + PI/2.0) * dPhi * dTheta;

    return integral;
}

double FullSteradian::get_integral() const
{
    double integral = 0;

    for (auto phi : this->spherical.phi)
        for (size_t theta_idx = 0; theta_idx < this->spherical.theta.size(); ++theta_idx)
            integral += sin(phi+PI/2.0) * dPhi * dTheta;

    return integral;
}

// Explicit instantiations for the types you want:
template double FullSteradian::get_integral<double>(std::vector<double>&) const;
template complex128 FullSteradian::get_integral<complex128>(std::vector<complex128>&) const;

template double FullSteradian::get_cos_integral<double>(const std::vector<double>&) const;
template complex128 FullSteradian::get_cos_integral<complex128>(const std::vector<complex128>&) const;
// --
