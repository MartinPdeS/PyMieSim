#pragma once


#include "utils.cpp"
#include "sources.cpp"
#include "fibonacci_mesh.cpp"
#include <cmath>
#include <vector>
#include <complex>
#include <base_class.cpp>

typedef std::complex<double> complex128;

class BaseScatterer : public Base
{
public:
    virtual ~BaseScatterer() = default;
    BaseScatterer() = default;

    BaseScatterer(const SOURCE::BaseSource &source, const double medium_index) : Base(max_order, source, medium_index){}

};

