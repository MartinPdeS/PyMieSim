#pragma once

#include "utils/base_mesh.h"

#define PI (double)3.14159265358979323846264338

class FullSteradian : public BaseMesh
{
    public:
        double dTheta, dPhi;

        FullSteradian(const size_t sampling, const double radius);

        template<typename dtype> dtype get_integral(std::vector<dtype>& vector) const;

        template<typename dtype> dtype get_cos_integral(const std::vector<dtype>& vector) const;

        double get_integral() const;

};

// --
