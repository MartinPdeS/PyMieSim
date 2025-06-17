#pragma once

#include "../coordinates/coordinates.h"


class BaseMesh
{
    public:
        size_t sampling;
        double radius;

    BaseMesh() = default;
    ~BaseMesh(){}

    BaseMesh(size_t sampling, double radius) : sampling(sampling), radius(radius) {}

};
