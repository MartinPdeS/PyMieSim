#include <cstddef>
#include <stdexcept>

// trapz over uniformly spaced samples: y[0..n-1], spacing dx
inline double trapz(std::vector<double> y, double dx) {
    const size_t n = y.size();
    if (n < 2)
        return 0.0;

    double s = 0.0;
    for (size_t i = 1; i < n; ++i)
        s += y[i - 1] + y[i];

    return 0.5 * dx * s;
}

// trapz over nonuniform x: sum_i 0.5*(x[i+1]-x[i])*(y[i]+y[i+1])
inline double trapz(std::vector<double> x, std::vector<double> y) {
    const size_t n = y.size();

    if (x.size() != n)
        throw std::invalid_argument("x and y must have same length");

    if (n < 2)
        return 0.0;

    double s = 0.0;
    for (size_t i = 0; i + 1 < n; ++i)
        s += 0.5 * (x[i + 1] - x[i]) * (y[i] + y[i + 1]);

    return s;
}
