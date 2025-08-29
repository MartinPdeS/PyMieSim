#include "coordinates.h"


// ---------- Vector Field ----------------
VectorField::VectorField(const std::vector<double>& vector)
: sampling(1), shape({1, 3}), data(vector) {}

VectorField::VectorField(size_t sampling)
: sampling(sampling), data(3 * sampling, 0.0) {}


// Accessors
double& VectorField::operator[](size_t i) { return data[i]; }
const double& VectorField::operator[](size_t i) const { return data[i]; }

double& VectorField::at(size_t i, size_t j) { return data[i * 3 + j]; }
const double& VectorField::at(size_t i, size_t j) const { return data[i * 3 + j]; }

double& VectorField::operator()(const size_t i, const size_t j) { return data[i * 3 + j]; }
const double& VectorField::operator()(const size_t i, const size_t j) const { return data[i * 3 + j]; }


// Vector addition with operator overloading
VectorField VectorField::operator+(const VectorField& other) const {
    VectorField result(sampling);
    for (size_t i = 0; i < data.size(); ++i) {
        result.data[i] = data[i] + other.data[i];
    }
    return result;
}

// Scalar product with a base vector
std::vector<double> VectorField::get_scalar_product(const VectorField& base) const {
    std::vector<double> output(sampling);
    for (size_t i = 0; i < sampling; ++i) {
        output[i] = (*this)(i, 0) * base[0] + (*this)(i, 1) * base[1] + (*this)(i, 2) * base[2];
    }
    return output;
}

// Normalize each vector in the field
void VectorField::normalize() {
    for (size_t i = 0; i < sampling; ++i) {
        double norm = std::sqrt(at(i, 0) * at(i, 0) + at(i, 1) * at(i, 1) + at(i, 2) * at(i, 2));
        if (norm > 0.0) {
            at(i, 0) /= norm;
            at(i, 1) /= norm;
            at(i, 2) /= norm;
        }
    }
}

// Rotate the vector field about a given axis by an angle
void VectorField::rotate_about_axis(const char axis, const double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    std::vector<std::vector<double>> matrix(3, std::vector<double>(3));

    switch (axis) {
        case 'x':
            matrix = {{1, 0, 0}, {0, c, -s}, {0, s, c}};
            break;
        case 'y':
            matrix = {{c, 0, s}, {0, 1, 0}, {-s, 0, c}};
            break;
        case 'z':
            matrix = {{c, -s, 0}, {s, c, 0}, {0, 0, 1}};
            break;
        default:
            throw std::invalid_argument("Invalid axis for rotation.");
    }

    apply_matrix(matrix);
}

// Apply a 3x3 rotation matrix to all vectors in the field
void VectorField::apply_matrix(const std::vector<std::vector<double>>& matrix) {
    std::array<double, 3> temp;
    for (size_t i = 0; i < sampling; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            temp[j] = matrix[j][0] * at(i, 0) + matrix[j][1] * at(i, 1) + matrix[j][2] * at(i, 2);
        }
        std::copy(temp.begin(), temp.end(), &data[i * 3]);
    }
}


// ---------- Spherical Coordinate System ----------------
Spherical::Spherical(const size_t sampling){
    r.reserve(sampling);
    phi.reserve(sampling);
    theta.reserve(sampling);
}

// ---------- Cartesian Coordinate System ----------------
Cartesian::Cartesian(const size_t sampling){
    x.reserve(sampling);
    y.reserve(sampling);
    z.reserve(sampling);
}

void Cartesian::set_coordinates_py(const std::vector<double> &x_values, const std::vector<double> &y_values, const std::vector<double> &z_values) {
    x = x_values;
    y = y_values;
    z = z_values;
}

Spherical Cartesian::to_spherical() const {
    Spherical sph(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        sph.r.push_back(std::sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]));
        sph.theta.push_back(std::atan2(y[i], x[i]));
        sph.phi.push_back(std::asin(z[i] / sph.r[i]));
    }
    return sph;
}

void Cartesian::rotate_about_axis(const char axis, const double angle) {
    const double c = std::cos(angle), s = std::sin(angle);
    std::vector<std::vector<double>> matrix;

    switch (axis) {
        case 'x':
            matrix = {{1, 0, 0}, {0, c, -s}, {0, s, c}};
            break;
        case 'y':
            matrix = {{c, 0, s}, {0, 1, 0}, {-s, 0, c}};
            break;
        case 'z':
            matrix = {{c, -s, 0}, {s, c, 0}, {0, 0, 1}};
            break;
        default:
            throw std::invalid_argument("Invalid axis for rotation.");
    }
    apply_matrix(matrix);
}

void Cartesian::apply_matrix(const std::vector<std::vector<double>>& matrix) {
    double tempx, tempy, tempz;
    for (size_t i = 0; i < x.size(); i++){

        tempx = matrix[0][0] * x[i] + matrix[0][1] * y[i] + matrix[0][2] * z[i];
        tempy = matrix[1][0] * x[i] + matrix[1][1] * y[i] + matrix[1][2] * z[i];
        tempz = matrix[2][0] * x[i] + matrix[2][1] * y[i] + matrix[2][2] * z[i];

        x[i] = tempx;
        y[i] = tempy;
        z[i] = tempz;
    }
}
