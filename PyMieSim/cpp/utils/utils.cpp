#pragma once

#include <algorithm>
#include <vector>
#include <complex>
#include <cmath>

typedef std::complex<double> complex128;
#define PI (double)3.14159265358979323846264338



// template <class T>
// std::vector <std::vector<T>> matrix_multiply(std::vector<std::vector<T>> &a, std::vector <std::vector<T>> &b)
// {
//     const int n = a.size();     // a rows
//     const int m = a[0].size();  // a cols
//     const int p = b[0].size();  // b cols

//     std::vector <std::vector<T>> c(n, std::vector<T>(p, 0));
//     for (auto j = 0; j < p; ++j)
//     {
//         for (auto k = 0; k < m; ++k)
//         {
//             for (auto i = 0; i < n; ++i)
//             {
//                 c[i][j] += a[i][k] * b[k][j];
//             }
//         }
//     }
//     return c;
// }

// std::vector<std::vector<double>> get_rotation_matrix(std::vector<double> rotation_axis, double rotation_angle)
// {
//     double norm_rotation_axis = sqrt(pow(rotation_axis[0], 2) + pow(rotation_axis[1], 2) + pow(rotation_axis[2], 2));

//     for (double &x: rotation_axis)
//         x /= norm_rotation_axis;

//     double
//         a = cos(rotation_angle / 2.0),
//         b = -1 * sin(rotation_angle / 2.0) * rotation_axis[0],
//         c = -1 * sin(rotation_angle / 2.0) * rotation_axis[1],
//         d = -1 * sin(rotation_angle / 2.0) * rotation_axis[2];

//     std::vector<std::vector<double>> matrix = {
//         {a * a + b * b - c * c - d * d, 2 * (b * c + a * d), 2 * (b * d - a * c)},
//         {2 * (b * c - a * d), a * a + c * c - b * b - d * d, 2 * (c * d + a * b)},
//         {2 * (b * d + a * c), 2 * (c * d - a * b), a * a + d * d - b * b - c * c}
//     };

//     return matrix;


// }
// -
