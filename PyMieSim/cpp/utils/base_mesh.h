#pragma once

#include "utils/coordinates.cpp"


class BaseMesh
{
    public:
        size_t sampling;
        double radius;
        Spherical spherical_coordinates;
        Cartesian cartesian_coordinates;

    BaseMesh() = default;
    ~BaseMesh(){}

    BaseMesh(size_t sampling, double radius)
    :   sampling(sampling),
        radius(radius),
        spherical_coordinates(sampling),
        cartesian_coordinates(sampling){}

};
