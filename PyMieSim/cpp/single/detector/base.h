#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <utils/math.h>
#include <utils/constants.h>
#include <material/material.h>
#include <mesh/fibonacci.h>

#include <single/source/source.h>
#include <single/optical_interface/optical_interface.h>
#include <single/scatterer/base_scatterer/base_scatterer.h>
#include <single/mode_field/mode_field.h>


using complex128 = std::complex<double>;


class BaseDetector {
public:
    struct CartesianCoordinateVectors {
        /**
         * @brief Cartesian coordinates of the detector angular sampling directions.
         *
         * The three vectors have length ``sampling`` and represent points on the
         * unit sphere associated with the detector Fibonacci mesh.
         */
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;
    };

    struct CollectionConeSurface {
        /**
         * @brief Cartesian surface coordinates for the detector collection cone.
         *
         * The three vectors are flattened row-major arrays with shape
         * ``radial_sampling`` by ``angular_sampling``. They are intended for
         * plotting utilities such as Matplotlib ``plot_surface``.
         */
        std::vector<double> x;
        std::vector<double> y;
        std::vector<double> z;

        std::size_t radial_sampling = 0;
        std::size_t angular_sampling = 0;
    };

    size_t sampling = 0;
    double numerical_aperture = 0.0;
    double cache_numerical_aperture = 0.0;
    double phi_offset = 0.0;
    double gamma_offset = 0.0;
    PolarizationState polarization_filter;
    std::shared_ptr<BaseMedium> medium;
    bool mean_coupling = true;

    std::vector<size_t> indices;
    double max_angle = 0.0;
    double min_angle = 0.0;
    std::vector<complex128> scalar_field;
    FibonacciMesh fibonacci_mesh;

    ModeID mode_id;
    ModeField mode_field;

    SnellFresnelInterface snell_interface;
    std::vector<double> interface_t_s;
    std::vector<double> interface_t_p;

    BaseDetector() = default;
    virtual ~BaseDetector() = default;

    BaseDetector(
        const size_t _sampling,
        const double _numerical_aperture,
        const double _cache_numerical_aperture,
        const double _phi_offset,
        const double _gamma_offset,
        const PolarizationState polarization_filter,
        std::shared_ptr<BaseMedium> _medium,
        const bool _mean_coupling
    )
    :   sampling(_sampling),
        numerical_aperture(_numerical_aperture),
        cache_numerical_aperture(_cache_numerical_aperture),
        phi_offset(_phi_offset),
        gamma_offset(_gamma_offset),
        polarization_filter(polarization_filter),
        medium(std::move(_medium)),
        mean_coupling(_mean_coupling)
    {}

    /**
     * @brief Initializes the detector angular mesh for a given scatterer.
     *
     * Implementations configure the Fibonacci mesh, detector angular limits,
     * scalar collecting field, and interface transmission coefficients required
     * for coupling calculations.
     *
     * @param scatterer Scatterer used to initialize detector-dependent optical
     *                  quantities.
     */
    virtual void initialize_mesh(const std::shared_ptr<BaseScatterer> scatterer) = 0;

    /**
     * @brief Computes the optical power coupled from a scatterer into the detector.
     *
     * @param scatterer Scatterer that produces the scattered field.
     * @param source Incident source illuminating the scatterer.
     *
     * @return Coupled optical power in watts before Python-side unit wrapping.
     */
    virtual double get_coupling(
        std::shared_ptr<BaseScatterer> scatterer,
        std::shared_ptr<BaseSource> source
    ) = 0;

    /**
     * @brief Parses a mode identifier string.
     *
     * The input is expected to follow the supported mode convention, for example
     * ``LP01`` or ``HG12``. The parsed result is stored in ``mode_id``.
     *
     * @param mode_number Mode identifier string.
     *
     * @throws std::invalid_argument If the mode identifier cannot be parsed.
     */
    virtual void parse_mode(const std::string& mode_number);

    /**
     * @brief Computes a structured representation of the detector scalar field.
     *
     * @param sampling Number of samples per structured angular dimension.
     *
     * @return Flattened complex scalar field with implementation-defined layout.
     */
    virtual std::vector<complex128> get_structured_scalarfield(const size_t sampling) const = 0;

    /**
     * @brief Returns the detector angular mesh as Cartesian coordinate vectors.
     *
     * The returned coordinates correspond to the current Fibonacci mesh and are
     * suitable for plotting the detector sampling directions on the unit sphere.
     *
     * @return Cartesian coordinate vectors ``x``, ``y``, and ``z``.
     */
    CartesianCoordinateVectors get_cartesian_coordinate_vectors() const {
        CartesianCoordinateVectors coordinates;

        coordinates.x.resize(this->sampling);
        coordinates.y.resize(this->sampling);
        coordinates.z.resize(this->sampling);

        for (std::size_t index = 0; index < this->sampling; ++index) {
            coordinates.x[index] = this->fibonacci_mesh.cartesian.x[index];
            coordinates.y[index] = this->fibonacci_mesh.cartesian.y[index];
            coordinates.z[index] = this->fibonacci_mesh.cartesian.z[index];
        }

        return coordinates;
    }

    /**
     * @brief Computes the mean collection axis of the detector mesh.
     *
     * The axis is obtained from the mean of the detector sampling directions and
     * normalized to unit length. It is primarily used to orient plotting geometry,
     * such as a collection cone.
     *
     * @return Unit vector representing the mean collection direction.
     *
     * @throws std::runtime_error If the detector mesh does not define a valid
     *                            mean direction.
     */
    std::array<double, 3> get_collection_axis() const {
        if (this->sampling == 0) {
            throw std::runtime_error("Cannot compute detector collection axis because sampling is zero.");
        }

        double axis_x = 0.0;
        double axis_y = 0.0;
        double axis_z = 0.0;

        for (std::size_t index = 0; index < this->sampling; ++index) {
            axis_x += this->fibonacci_mesh.cartesian.x[index];
            axis_y += this->fibonacci_mesh.cartesian.y[index];
            axis_z += this->fibonacci_mesh.cartesian.z[index];
        }

        const double axis_norm = std::sqrt(
            axis_x * axis_x +
            axis_y * axis_y +
            axis_z * axis_z
        );

        if (axis_norm <= 0.0) {
            throw std::runtime_error("Cannot compute detector collection axis because the mean mesh direction is zero.");
        }

        return {
            axis_x / axis_norm,
            axis_y / axis_norm,
            axis_z / axis_norm
        };
    }

    /**
     * @brief Builds a cone surface representing the detector collection aperture.
     *
     * The cone is oriented along the detector collection axis and opens with
     * half-angle ``max_angle``. The returned surface is flattened in row-major
     * order with shape ``radial_sampling`` by ``angular_sampling``.
     *
     * This method only generates geometry. Rendering is intentionally left to
     * the Python binding layer so that plotting backends remain decoupled from
     * the detector model.
     *
     * @param angular_sampling Number of angular samples used around the cone.
     * @param radial_sampling Number of radial samples used from the apex to the rim.
     *
     * @return Cartesian surface coordinates for the collection cone.
     *
     * @throws std::invalid_argument If either sampling argument is zero.
     */
    CollectionConeSurface get_collection_cone_surface(
        const std::size_t angular_sampling = 96,
        const std::size_t radial_sampling = 16
    ) const {
        if (angular_sampling == 0) {
            throw std::invalid_argument("Cone angular sampling must be greater than zero.");
        }

        if (radial_sampling == 0) {
            throw std::invalid_argument("Cone radial sampling must be greater than zero.");
        }

        const double pi = 3.14159265358979323846;

        const std::array<double, 3> axis = this->get_collection_axis();

        double reference_x = 0.0;
        double reference_y = 0.0;
        double reference_z = 1.0;

        if (std::abs(axis[2]) > 0.9) {
            reference_x = 1.0;
            reference_y = 0.0;
            reference_z = 0.0;
        }

        double basis_x_x = reference_y * axis[2] - reference_z * axis[1];
        double basis_x_y = reference_z * axis[0] - reference_x * axis[2];
        double basis_x_z = reference_x * axis[1] - reference_y * axis[0];

        const double basis_x_norm = std::sqrt(
            basis_x_x * basis_x_x +
            basis_x_y * basis_x_y +
            basis_x_z * basis_x_z
        );

        if (basis_x_norm <= 0.0) {
            throw std::runtime_error("Cannot build detector cone because the transverse basis is degenerate.");
        }

        basis_x_x /= basis_x_norm;
        basis_x_y /= basis_x_norm;
        basis_x_z /= basis_x_norm;

        const double basis_y_x = axis[1] * basis_x_z - axis[2] * basis_x_y;
        const double basis_y_y = axis[2] * basis_x_x - axis[0] * basis_x_z;
        const double basis_y_z = axis[0] * basis_x_y - axis[1] * basis_x_x;

        CollectionConeSurface surface;
        surface.angular_sampling = angular_sampling;
        surface.radial_sampling = radial_sampling;

        const std::size_t surface_size = angular_sampling * radial_sampling;

        surface.x.resize(surface_size);
        surface.y.resize(surface_size);
        surface.z.resize(surface_size);

        const double cos_max_angle = std::cos(this->max_angle);
        const double sin_max_angle = std::sin(this->max_angle);

        for (std::size_t radial_index = 0; radial_index < radial_sampling; ++radial_index) {
            const double radius_fraction = static_cast<double>(radial_index) /
                static_cast<double>(std::max<std::size_t>(radial_sampling - 1, 1));

            for (std::size_t angular_index = 0; angular_index < angular_sampling; ++angular_index) {
                const double theta = 2.0 * pi *
                    static_cast<double>(angular_index) /
                    static_cast<double>(std::max<std::size_t>(angular_sampling - 1, 1));

                const double rim_x =
                    cos_max_angle * axis[0] +
                    sin_max_angle * (
                        std::cos(theta) * basis_x_x +
                        std::sin(theta) * basis_y_x
                    );

                const double rim_y =
                    cos_max_angle * axis[1] +
                    sin_max_angle * (
                        std::cos(theta) * basis_x_y +
                        std::sin(theta) * basis_y_y
                    );

                const double rim_z =
                    cos_max_angle * axis[2] +
                    sin_max_angle * (
                        std::cos(theta) * basis_x_z +
                        std::sin(theta) * basis_y_z
                    );

                const std::size_t flat_index = radial_index * angular_sampling + angular_index;

                surface.x[flat_index] = radius_fraction * rim_x;
                surface.y[flat_index] = radius_fraction * rim_y;
                surface.z[flat_index] = radius_fraction * rim_z;
            }
        }

        return surface;
    }

    /**
     * @brief Squares each element in the input vector in place.
     *
     * @tparam T Element type.
     * @param array Input vector to modify.
     */
    template <typename T> inline void square_array(std::vector<T>& array);

    /**
     * @brief Computes the squared absolute value of the sum of vector elements.
     *
     * @tparam T Element type.
     * @param array Input vector.
     *
     * @return Squared one-norm-like coherent sum.
     */
    template <typename T> inline double get_norm1_squared(const std::vector<T>& array) const;

    /**
     * @brief Applies cached interface transmission coefficients to angular fields.
     *
     * @param theta_field Field component along the local theta direction.
     * @param phi_field Field component along the local phi direction.
     */
    void apply_interface_transmission_to_fields(
        std::vector<complex128>& theta_field,
        std::vector<complex128>& phi_field
    ) const;

    /**
     * @brief Computes the squared two-norm of a vector.
     *
     * @tparam T Element type.
     * @param array Input vector.
     *
     * @return Sum of squared absolute values.
     */
    template <typename T> inline double get_norm2_squared(const std::vector<T>& array) const;

    /**
     * @brief Computes Poynting field samples for a scatterer and source.
     *
     * @param scatterer Scatterer that produces the scattered field.
     * @param source Incident source illuminating the scatterer.
     * @param distance Observation distance in meters.
     *
     * @return Poynting vector magnitudes sampled on the detector mesh.
     */
    std::vector<double> get_poynting_field(
        std::shared_ptr<BaseScatterer> scatterer,
        std::shared_ptr<BaseSource> source,
        double distance = 1.0
    ) const;

    /**
     * @brief Computes the total radiated energy flow sampled by the detector mesh.
     *
     * @param scatterer Scatterer that produces the scattered field.
     * @param source Incident source illuminating the scatterer.
     *
     * @return Total energy flow in watts before Python-side unit wrapping.
     */
    double get_energy_flow(
        std::shared_ptr<BaseScatterer> scatterer,
        std::shared_ptr<BaseSource> source
    ) const;

    /**
     * @brief Prints detector properties to standard output.
     *
     * @param precision Number of decimal places used for floating point values.
     */
    void print_properties(int precision = 4) const {
        std::printf("Detector Properties:\n");
        std::printf("  Sampling: %zu\n", this->sampling);
        std::printf("  Numerical Aperture: %.*f\n", precision, this->numerical_aperture);
        std::printf("  Cache Numerical Aperture: %.*f\n", precision, this->cache_numerical_aperture);
        std::printf("  Phi Offset (radians): %.*f\n", precision, this->phi_offset);
        std::printf("  Gamma Offset (radians): %.*f\n", precision, this->gamma_offset);
        std::printf("  Polarization Filter: %.*f\n", precision, this->polarization_filter.angle);
    }

private:
    static inline double clamp_m1_p1(double value) {
        return std::max(-1.0, std::min(1.0, value));
    }
};