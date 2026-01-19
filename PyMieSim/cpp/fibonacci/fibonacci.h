#pragma once


#include <vector>
#include <cmath>
#include "../coordinates/coordinates.h"
#include <utils/constants.h>


class FibonacciMesh {
    public:
        size_t sampling;
        size_t true_number_of_sample;
        double max_angle = 0.0;
        double min_angle = 0.0;
        double phi_offset = 0.0;
        double gamma_offset = 0.0;
        double dOmega = 0.0;
        double Omega = 0.0;
        double radius;
        double rotation;

        Cartesian base_cartesian;
        Cartesian cartesian;
        Spherical spherical;

        VectorField perpendicular_vector;
        VectorField parallel_vector;
        VectorField vertical_vector_field{{1, 0, 0}};
        VectorField horizontal_vector_field{{0, 1, 0}};

        std::vector<double> horizontal_parallel_projection;
        std::vector<double> vertical_parallel_projection;
        std::vector<double> horizontal_perpendicular_projection;
        std::vector<double> vertical_perpendicular_projection;

        FibonacciMesh() = default;

        /**
         *  @brief Constructs a Fibonacci mesh with the specified parameters.
         *  @param sampling The number of sampling points in the mesh.
         *  @param max_angle The maximum angle for the mesh points (in radians).
         *  @param min_angle The minimum angle for the mesh points (in radians).
         *  @param phi_offset The offset for the phi angle (in radians).
         *  @param gamma_offset The offset for the gamma angle (in radians).
         *  @param rotation The rotation angle for the mesh (in radians).
         *  @param radius The radius of the mesh (default is 1.0).
         */
        FibonacciMesh(size_t sampling, double max_angle, double min_angle, double phi_offset, double gamma_offset, double rotation, double radius = 1.0);

        /**
         *  @brief Rotates the mesh around its center based on the specified phi and gamma offsets.
         *  @note This function applies rotations to the Cartesian coordinates and vector fields.
         */
        void rotate_around_center();

        /**
         *  @brief Computes the vector field for the Fibonacci mesh.
         *  @note This function initializes the parallel and perpendicular vector fields based on the horizontal and vertical vector fields.
         */
        void compute_vector_field();

        /**
         *  @brief Computes the projections of the vector fields onto the parallel and perpendicular vectors.
         *  @note This function calculates the projections based on the horizontal and vertical vector fields.
         */
        void compute_projections();

        /**
         *  @brief Rotates the mesh around the specified axis by the given angle.
         *  @param angle The angle by which to rotate (in radians).
         */
        void rotate_around_axis(double angle);

        /**
         *  @brief Computes the Fibonacci mesh based on the specified parameters.
         *  @note This function generates the mesh points in Cartesian coordinates and computes the corresponding spherical coordinates.
         */
        void compute_mesh();

        /**
         *  @brief Computes the properties of the Fibonacci mesh, including solid angle and sampling ratio.
         */
        void compute_properties();

        /**
         *  @brief Gets the principal axis of the Fibonacci mesh.
         *  @return A vector containing the principal axis coordinates.
         */
        std::vector<double> get_principal_axis() const;

        /**
         *  @brief Computes the rotation matrix for a given rotation axis and angle.
         *  @param rotation_axis The axis around which to rotate.
         *  @param rotation_angle The angle by which to rotate (in radians).
         *  @return The rotation matrix as a 3x3 vector.
         */
        std::vector<std::vector<double>> get_rotation_matrix(std::vector<double> rotation_axis, double rotation_angle) const;

};
