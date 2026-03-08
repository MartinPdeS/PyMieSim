#pragma once

#include <vector>
#include <complex>
#include <string>
#include <stdexcept>

typedef std::complex<double> complex128;

class DispersiveMaterial {
    std::vector<complex128> refractive_indices;

    DispersiveMaterial = default;
    DispersiveMaterial(const std::vector<complex128>& refractive_indices) : refractive_indices(refractive_indices) {}

};