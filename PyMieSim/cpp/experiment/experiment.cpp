#include <cstdarg>    // For va_list, va_start, va_end
#include <cstdio>     // For printf, vprintf
#include "sphere/sphere.h"
#include "cylinder/cylinder.h"
#include "coreshell/coreshell.h"
#include "source/source.h"
#include "detector/detector.h"
#include "utils/numpy_interface.h"
#include "sets/sets.cpp"


#define DEFINE_SCATTERER_FUNCTION(scatterer, dtype, name) \
    pybind11::array_t<dtype> get_##scatterer##_##name() const { \
        debug_printf("Computing %s_%s\n", #scatterer, #name); \
        return get_scatterer_data<double, scatterer##Set>(scatterer##_set, &scatterer::get_##name); \
    } \
    pybind11::array_t<dtype> get_##scatterer##_##name##_sequential() const { \
        debug_printf("Computing %s_%s_sequential\n", #scatterer, #name); \
        return get_scatterer_data_sequential<double, scatterer##Set>(scatterer##_set, &scatterer::get_##name); \
    }

#define DEFINE_SCATTERER_COEFFICIENT(scatterer, name) \
    pybind11::array_t<double> get_##scatterer##_##name() const { \
        debug_printf("Computing %s_%s coefficient\n", #scatterer, #name); \
        return get_scatterer_data<double, scatterer##Set>(scatterer##_set, &scatterer::get_##name##_abs); \
    } \
    pybind11::array_t<double> get_##scatterer##_##name##_sequential() const { \
        debug_printf("Computing %s_%s coefficient sequentially\n", #scatterer, #name); \
        return get_scatterer_data_sequential<double, scatterer##Set>(scatterer##_set, &scatterer::get_##name##_abs); \
    }

#define DEFINE_SCATTERER_COUPLING(scatterer) \
    pybind11::array_t<double> get_##scatterer##_coupling() const { \
        debug_printf("Computing %s_coupling\n", #scatterer); \
        return get_scatterer_coupling<scatterer##Set>(scatterer##_set); \
    } \
    pybind11::array_t<double> get_##scatterer##_coupling_sequential() const { \
        debug_printf("Computing %s_coupling_sequential\n", #scatterer); \
        return get_scatterer_coupling_sequential<scatterer##Set>(scatterer##_set); \
    }

class Experiment
{
    public:
        bool debug_mode = true;
        SphereSet Sphere_set;
        CylinderSet Cylinder_set;
        CoreShellSet CoreShell_set;
        DetectorSet Detector_set;
        SourceSet Source_set;

        explicit Experiment(bool debug_mode = true) : debug_mode(debug_mode) {}

        // Setter methods
        void set_sphere(SphereSet& set) { Sphere_set = set; }
        void set_cylinder(CylinderSet& set) { Cylinder_set = set; }
        void set_coreshell(CoreShellSet& set) { CoreShell_set = set; }
        void set_source(SourceSet &set) { Source_set = set; }
        void set_detector(DetectorSet &set) { Detector_set = set; }

        // Helper method for debugging without iostream.
        // It prints the given message only if debug_mode is true.
        void debug_printf(const char* format, ...) const {
            if (!debug_mode) return;
            va_list args;
            va_start(args, format);
            vprintf(format, args);
            va_end(args);
        }

        // Example: You can also add a class-wide debug_print method.
        void debug_print_state() const {
            if (!debug_mode) return;
            debug_printf("----- Experiment Debug Info -----\n");
            debug_printf("SourceSet total combinations: %zu\n", Source_set.total_combinations);
            debug_printf("SphereSet total combinations: %zu\n", Sphere_set.total_combinations);
            debug_printf("CylinderSet total combinations: %zu\n", Cylinder_set.total_combinations);
            debug_printf("CoreshellSet total combinations: %zu\n", CoreShell_set.total_combinations);
            debug_printf("DetectorSet total combinations: %zu\n", Detector_set.total_combinations);
            debug_printf("SourceSet shape: ");
            for (size_t dim : Source_set.shape) {
                debug_printf("%zu ", dim);
            }
            debug_printf("\n---------------------------------\n");
        }

        // Flattens a multi-index into a single index.
        template <typename... MultiIndexVectors>
        size_t flatten_multi_index(const std::vector<size_t>& dimensions, const MultiIndexVectors&... multi_indices) const {
            std::vector<size_t> multi_index = concatenate_vector(multi_indices...);
            size_t flatten_index = 0;
            size_t stride = 1;
            for (int i = static_cast<int>(dimensions.size()) - 1; i >= 0; --i) {
                flatten_index += multi_index[i] * stride;
                stride *= dimensions[i];
            }
            return flatten_index;
        }

        template <typename T, typename... Ts>
        T concatenate_vector(const T& first_vector, const Ts&... other_vectors) const
        {
            T output_vector = first_vector;
            (output_vector.insert(output_vector.end(), other_vectors.begin(), other_vectors.end()), ...);
            return output_vector;
        }

        template<typename dtype, typename ScattererSet, typename Function>
        pybind11::array_t<dtype> get_scatterer_data(const ScattererSet& scatterer_set, Function function) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape;
            size_t total_iterations;

            if (Detector_set.is_empty) {
                array_shape = concatenate_vector(Source_set.shape, scatterer_set.shape);
                total_iterations = Source_set.total_combinations * scatterer_set.total_combinations;
            }
            else {
                array_shape = concatenate_vector(Source_set.shape, scatterer_set.shape, Detector_set.shape);
                total_iterations = Source_set.total_combinations * scatterer_set.total_combinations * Detector_set.total_combinations;
            }

            debug_printf("get_scatterer_data: total_iterations = %zu\n", total_iterations);

            std::vector<dtype> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t flat_index = 0; flat_index < total_iterations; ++flat_index) {
                size_t idx; // Declare idx locally so each iteration has its own copy

                if (Detector_set.is_empty) {
                    // 2D case: only source and scatterer
                    size_t i = flat_index / scatterer_set.total_combinations;
                    size_t j = flat_index % scatterer_set.total_combinations;
                    BaseSource source = Source_set.get_source_by_index(i);
                    auto scatterer = scatterer_set.get_scatterer_by_index(j, source);
                    idx = flatten_multi_index(array_shape, source.indices, scatterer.indices);
                    dtype value = (scatterer.*function)();
                    output_array[idx] = value;
                } else {
                    // 3D case: source, scatterer, and detector
                    size_t i = flat_index / (scatterer_set.total_combinations * Detector_set.total_combinations);
                    size_t j = (flat_index / Detector_set.total_combinations) % scatterer_set.total_combinations;
                    size_t k = flat_index % Detector_set.total_combinations;
                    BaseSource source = Source_set.get_source_by_index(i);
                    auto scatterer = scatterer_set.get_scatterer_by_index(j, source);
                    Detector detector = Detector_set.get_detector_by_index(k);
                    idx = flatten_multi_index(array_shape, source.indices, scatterer.indices, detector.indices);
                    dtype value = (scatterer.*function)();
                    output_array[idx] = value;
                }
            }

            debug_printf("get_scatterer_data: finished computation\n");
            return _vector_to_numpy(output_array, array_shape);
        }

        template<typename dtype, typename ScattererSet, typename Function>
        pybind11::array_t<dtype> get_scatterer_data_sequential(const ScattererSet& scatterer_set, Function function) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape = {Source_set.wavelength.size()};
            size_t full_size = Source_set.wavelength.size();
            scatterer_set.validate_sequential_data(full_size);
            Source_set.validate_sequential_data(full_size);
            debug_printf("get_scatterer_data_sequential: full_size = %zu\n", full_size);

            std::vector<dtype> output_array(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {
                BaseSource source = Source_set.get_source_by_index_sequential(idx);
                auto scatterer = scatterer_set.get_scatterer_by_index_sequential(idx, source);
                output_array[idx] = (scatterer.*function)();
            }
            debug_printf("get_scatterer_data_sequential: finished computation\n");
            return _vector_to_numpy(output_array, {full_size});
        }

        template<typename ScattererSet>
        pybind11::array_t<double> get_scatterer_coupling(const ScattererSet& scatterer_set) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape = concatenate_vector(Source_set.shape, scatterer_set.shape, Detector_set.shape);
            size_t total_iterations = Source_set.total_combinations * scatterer_set.total_combinations * Detector_set.total_combinations;
            debug_printf("get_scatterer_coupling: total_iterations = %zu\n", total_iterations);

            std::vector<double> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t idx_flat = 0; idx_flat < total_iterations; ++idx_flat) {
                size_t i = idx_flat / (scatterer_set.total_combinations * Detector_set.total_combinations);
                size_t j = (idx_flat / Detector_set.total_combinations) % scatterer_set.total_combinations;
                size_t k = idx_flat % Detector_set.total_combinations;
                BaseSource source = Source_set.get_source_by_index(i);
                auto scatterer = scatterer_set.get_scatterer_by_index(j, source);
                Detector detector = Detector_set.get_detector_by_index(k);
                detector.medium_refractive_index = scatterer.medium_refractive_index;
                size_t idx = flatten_multi_index(array_shape, source.indices, scatterer.indices, detector.indices);
                output_array[idx] = detector.get_coupling(scatterer);
            }
            debug_printf("get_scatterer_coupling: finished computation\n");
            return _vector_to_numpy(output_array, array_shape);
        }

        template<typename ScattererSet>
        pybind11::array_t<double> get_scatterer_coupling_sequential(const ScattererSet& scatterer_set) const {

            if (debug_mode)
                this->debug_print_state();

            std::vector<size_t> array_shape = {Source_set.wavelength.size()};
            size_t full_size = Source_set.wavelength.size();
            scatterer_set.validate_sequential_data(full_size);
            Source_set.validate_sequential_data(full_size);
            Detector_set.validate_sequential_data(full_size);
            debug_printf("get_scatterer_coupling_sequential: full_size = %zu\n", full_size);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {
                BaseSource source = Source_set.get_source_by_index_sequential(idx);
                auto scatterer = scatterer_set.get_scatterer_by_index_sequential(idx, source);
                Detector detector = Detector_set.get_detector_by_index_sequential(idx);
                output_array[idx] = detector.get_coupling(scatterer);
            }
            debug_printf("get_scatterer_coupling_sequential: finished computation\n");
            return _vector_to_numpy(output_array, {full_size});
        }

        //--------------------------------------SPHERE------------------------------------
        DEFINE_SCATTERER_COEFFICIENT(Sphere, a1)
        DEFINE_SCATTERER_COEFFICIENT(Sphere, a2)
        DEFINE_SCATTERER_COEFFICIENT(Sphere, a3)
        DEFINE_SCATTERER_COEFFICIENT(Sphere, b1)
        DEFINE_SCATTERER_COEFFICIENT(Sphere, b2)
        DEFINE_SCATTERER_COEFFICIENT(Sphere, b3)

        DEFINE_SCATTERER_FUNCTION(Sphere, double, Qsca)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Qext)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Qabs)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Qpr)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Qback)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Qforward)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Qratio)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Csca)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Cext)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Cabs)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Cpr)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Cback)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Cratio)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, Cforward)
        DEFINE_SCATTERER_FUNCTION(Sphere, double, g)

        DEFINE_SCATTERER_COUPLING(Sphere)

        //--------------------------------------CYLINDER------------------------------------
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, a11)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, a12)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, a13)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, a21)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, a22)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, a23)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, b11)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, b12)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, b13)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, b21)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, b22)
        DEFINE_SCATTERER_COEFFICIENT(Cylinder, b23)

        DEFINE_SCATTERER_FUNCTION(Cylinder, double, Qsca)
        DEFINE_SCATTERER_FUNCTION(Cylinder, double, Qext)
        DEFINE_SCATTERER_FUNCTION(Cylinder, double, Qabs)
        DEFINE_SCATTERER_FUNCTION(Cylinder, double, Csca)
        DEFINE_SCATTERER_FUNCTION(Cylinder, double, Cext)
        DEFINE_SCATTERER_FUNCTION(Cylinder, double, Cabs)
        DEFINE_SCATTERER_FUNCTION(Cylinder, double, g)

        DEFINE_SCATTERER_COUPLING(Cylinder)

        //--------------------------------------CORESHELL------------------------------------
        DEFINE_SCATTERER_COEFFICIENT(CoreShell, a1)
        DEFINE_SCATTERER_COEFFICIENT(CoreShell, a2)
        DEFINE_SCATTERER_COEFFICIENT(CoreShell, a3)
        DEFINE_SCATTERER_COEFFICIENT(CoreShell, b1)
        DEFINE_SCATTERER_COEFFICIENT(CoreShell, b2)
        DEFINE_SCATTERER_COEFFICIENT(CoreShell, b3)

        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Qsca)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Qext)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Qabs)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Qpr)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Qback)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Qforward)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Qratio)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Csca)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Cext)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Cabs)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Cpr)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Cback)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Cratio)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, Cforward)
        DEFINE_SCATTERER_FUNCTION(CoreShell, double, g)

        DEFINE_SCATTERER_COUPLING(CoreShell)
};
