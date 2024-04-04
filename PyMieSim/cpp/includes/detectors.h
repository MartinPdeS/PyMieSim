#include "fibonacci_mesh.cpp"
#include "utils.cpp"

#include <vector>
#include <complex>
#include <cmath> // For std::isnan and std::pow

namespace DETECTOR {

    using complex128 = std::complex<double>;

    struct State {
        std::vector<complex128> scalar_field;
        double NA = 0.0, phi_offset = 0.0, gamma_offset = 0.0, polarization_filter = 0.0, rotation_angle = 0.0;
        bool coherent = true, point_coupling = true;

        State() = default;

        State(const std::vector<complex128>& scalar_field, double NA, double phi_offset, double gamma_offset,
              double polarization_filter, double rotation_angle, bool coherent, bool point_coupling)
            : scalar_field(scalar_field), NA(NA), phi_offset(phi_offset), gamma_offset(gamma_offset),
              polarization_filter(polarization_filter), rotation_angle(rotation_angle),
              coherent(coherent), point_coupling(point_coupling) {}
    };

    class Detector {
    public:
        FibonacciMesh fibonacci_mesh;
        State state;
        bool coherent = true;

        Detector() = default;

        Detector(const std::vector<complex128>& scalar_field, double NA, double phi_offset,
                 double gamma_offset, double polarization_filter, double rotation_angle, bool coherent, bool point_coupling)
            : state(scalar_field, NA, phi_offset, gamma_offset, polarization_filter, rotation_angle, coherent, point_coupling) {
            fibonacci_mesh = FibonacciMesh(
                scalar_field.size(), NA2Angle(state.NA), state.phi_offset, state.gamma_offset, state.rotation_angle);
        }

        Detector(const State& state) : state(state) {
            fibonacci_mesh = FibonacciMesh(
                state.scalar_field.size(), NA2Angle(state.NA), state.phi_offset, state.gamma_offset, state.rotation_angle);
        }

        void rotate_around_axis(double rotation_angle) {
            fibonacci_mesh.rotate_around_axis(rotation_angle);
        }

        // Template method for coupling calculations
        template <typename T>
        double get_coupling(T& scatterer) {
            // Method logic simplified with conditional operator
            if (state.coherent) {
                return state.point_coupling ? get_coupling_point_coherent(scatterer) : get_coupling_mean_coherent(scatterer);
            } else {
                return state.point_coupling ? get_coupling_point_no_coherent(scatterer) : get_coupling_mean_no_coherent(scatterer);
            }
        }

        // Implementation details for coupling calculations should be defined here
        template <typename T> double get_coupling_point_no_coherent(T& scatterer);
        template <typename T> double get_coupling_mean_no_coherent(T& scatterer);
        template <typename T> double get_coupling_point_coherent(T& scatterer);
        template <typename T> double get_coupling_mean_coherent(T& scatterer);

    private:
        // Private helper methods for coupling calculations
        template <typename T> double calculate_coupling(const T& scatterer, bool point, bool coherent);
        std::tuple<std::vector<complex128>, std::vector<complex128>> get_projected_fields(const std::vector<complex128>& theta_field, const std::vector<complex128>& phi_field) const;
        void apply_scalar_field(std::vector<complex128>& field0, std::vector<complex128>& field1) const;
        template <typename T> inline double get_norm1_squared(const std::vector<T>& array) const;
        template <typename T> inline double get_norm2_squared(const std::vector<T>& array) const;
        template <typename T> inline void square_array(std::vector<T>& array);
        template <typename T> inline void apply_polarization_filter(T& coupling_theta, T& coupling_phi, double polarization_filter) const;
    };


} // namespace DETECTOR


