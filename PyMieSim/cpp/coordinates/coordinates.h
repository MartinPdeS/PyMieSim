#pragma once
#include <vector>
#include <cmath>
#include <array>
#include <stdexcept>


struct SphericalCoordinate {
    double r = 0.0, phi = 0.0, theta = 0.0;
};

struct VectorField {
    size_t sampling = 0;
    std::vector<size_t> shape{1, 3};
    std::vector<double> data;

    /**
     * @brief Default constructor for VectorField.
     * Initializes an empty vector field.
     */
    VectorField() = default;

    /**
     * @brief Constructor that initializes the vector field with a given vector.
     * @param vector A vector containing the initial values for the vector field.
     */
    VectorField(const std::vector<double>& vector);

    /**
     * @brief Constructor that initializes the vector field with a specified sampling size.
     * @param sampling The number of samples to reserve for the vector field.
     */
    VectorField(size_t sampling);

    /**
     * @brief Access the vector field data using an index.
     * @param i Index of the vector component.
     */
    double& operator[](size_t i);

    /**
     * @brief Access the vector field data using an index.
     * @param i Index of the vector component.
     */
    const double& operator[](size_t i) const;

    /**
     * @brief Access the vector field data at a specific index and component.
     * @param i Index of the vector.
     * @param j Component index (0, 1, or 2).
     */
    double& at(size_t i, size_t j);

    /**
     * @brief Access the vector field data at a specific index and component.
     * @param i Index of the vector.
     * @param j Component index (0, 1, or 2).
     */
    const double& at(size_t i, size_t j) const;

    /**
     * @brief Access the vector field data using a tuple of indices.
     * @param i Index of the vector.
     * @param j Component index (0, 1, or 2).
     */
    double& operator()(const size_t i, const size_t j);

    /**
     * @brief Access the vector field data using a tuple of indices.
     * @param i Index of the vector.
     * @param j Component index (0, 1, or 2).
     */
    const double& operator()(const size_t i, const size_t j) const;

    /**
     * @brief Add another vector field to this one.
     * @param other The vector field to add.
     * @return A new VectorField that is the sum of this and the other.
     */
    VectorField operator+(const VectorField& other) const;

    /**
     * @brief Compute the scalar product with a base vector.
     * @param base The base vector to compute the scalar product with.
     * @return A vector containing the scalar products.
     */
    std::vector<double> get_scalar_product(const VectorField& base) const;

    /**
     * @brief Normalize each vector in the field.
     * @note This modifies the vectors in place, ensuring each vector has a unit length.
     */
    void normalize();

    /**
     * @brief Rotate the vector field about a given axis by an angle.
     * @param axis The axis to rotate about ('x', 'y', or 'z').
     * @param angle The angle to rotate by (in radians).
     */
    void rotate_about_axis(const char axis, const double angle);

    /**
     * @brief Apply a 3x3 rotation matrix to all vectors in the field.
     * @param matrix The rotation matrix to apply.
     */
    void apply_matrix(const std::vector<std::vector<double>>& matrix);
};

struct Spherical {
    std::vector<double> r, phi, theta;

    /**
     * @brief Default constructor for Spherical coordinates.
     * Initializes empty vectors for r, phi, and theta coordinates.
     */
    Spherical() = default;

    /**
     * @brief Constructor that initializes the Spherical coordinates with a specified sampling size.
     * @param sampling The number of samples to reserve for the Spherical coordinates.
     */
    explicit Spherical(const size_t sampling);
};

struct Cartesian {
    std::vector<double> x, y, z;

    /**
     * @brief Default constructor for Cartesian coordinates.
     * Initializes empty vectors for x, y, and z coordinates.
     */
    Cartesian() = default;

    /**
     * @brief Constructor that initializes the Cartesian coordinates with a specified sampling size.
     * @param sampling The number of samples to reserve for the Cartesian coordinates.
     */
    explicit Cartesian(const size_t sampling);

    /**
     * @brief Set the Cartesian coordinates from Python.
     * @param x_values Vector of x coordinates.
     * @param y_values Vector of y coordinates.
     * @param z_values Vector of z coordinates.
     */
    void set_coordinates_py(const std::vector<double> &x_values, const std::vector<double> &y_values, const std::vector<double> &z_values);

    /**
     * @brief Convert Cartesian coordinates to spherical coordinates.
     * @return A Spherical object containing the converted coordinates.
     */
    Spherical to_spherical() const;

    /**
     * @brief Rotate the Cartesian coordinates about a specified axis by a given angle.
     * @param axis The axis to rotate about ('x', 'y', or 'z').
     * @param angle The angle to rotate by (in radians).
     */
    void rotate_about_axis(const char axis, const double angle);

    /**
     * @brief Apply a 3x3 rotation matrix to all vectors in the field.
     * @param matrix The rotation matrix to apply.
     */
    void apply_matrix(const std::vector<std::vector<double>>& matrix);

};
