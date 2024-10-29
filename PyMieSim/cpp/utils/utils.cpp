#pragma once

#include <algorithm>
#include <vector>
#include <complex>
#include <cmath>

typedef std::complex<double> complex128;
#define PI (double)3.14159265358979323846264338

double NA2Angle(const double &NA)
{
    if (NA <= 1.0)
        return asin(NA);

    if (NA >= 1.0)
        return asin(NA-1.0) + PI/2.0;

    return 1.0;
}

template <typename T, typename... Ts>
T concatenate_vector(const T& first_vector, const Ts&... other_vectors)
{
    T output_vector = first_vector;
    (output_vector.insert(output_vector.end(), other_vectors.begin(), other_vectors.end()), ...);
    return output_vector;
}

template <typename T>
T get_vector_sigma(const std::vector<T> &vector)
{
    T sigma = 1;
    for (auto e: vector)
      sigma *= e;

    return sigma;
}

template <class T>
void Squared(std::vector<T>& vector)
{
    for (size_t iter=0; iter<vector.size(); iter++)
        vector[iter] = pow(abs(vector[iter]), 2);
}

template <class T>
std::vector<T> Add(std::vector<T>& vector0, std::vector<T>& vector1)
{
    std::vector<T> output_vector;
    output_vector.reserve( vector0.size() );

    for (size_t iter=0; iter<vector0.size(); iter++)
        output_vector.push_back( vector0[iter] + vector1[iter] );

    return output_vector;
}

template <class T>
std::vector <std::vector<T>> matrix_multiply(std::vector<std::vector<T>> &a, std::vector <std::vector<T>> &b)
{
    const int n = a.size();     // a rows
    const int m = a[0].size();  // a cols
    const int p = b[0].size();  // b cols

    std::vector <std::vector<T>> c(n, std::vector<T>(p, 0));
    for (auto j = 0; j < p; ++j)
    {
        for (auto k = 0; k < m; ++k)
        {
            for (auto i = 0; i < n; ++i)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return c;
}

std::vector<std::vector<double>> get_rotation_matrix(std::vector<double> rotation_axis, double rotation_angle)
{
    double norm_rotation_axis = sqrt(pow(rotation_axis[0], 2) + pow(rotation_axis[1], 2) + pow(rotation_axis[2], 2));

    for (double &x: rotation_axis)
        x /= norm_rotation_axis;

    double
        a = cos(rotation_angle / 2.0),
        b = -1 * sin(rotation_angle / 2.0) * rotation_axis[0],
        c = -1 * sin(rotation_angle / 2.0) * rotation_axis[1],
        d = -1 * sin(rotation_angle / 2.0) * rotation_axis[2];

    std::vector<std::vector<double>> matrix = {
        {a * a + b * b - c * c - d * d, 2 * (b * c + a * d), 2 * (b * d - a * c)},
        {2 * (b * c - a * d), a * a + c * c - b * b - d * d, 2 * (c * d + a * b)},
        {2 * (b * d + a * c), 2 * (c * d - a * b), a * a + d * d - b * b - c * c}
    };

    return matrix;


}
// -
