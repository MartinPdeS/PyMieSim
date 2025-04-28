#include <cstdarg>    // For va_list, va_start, va_end
#include <cstdio>     // For printf, vprintf
#include "scatterer/sphere/sphere.h"
#include "scatterer/cylinder/cylinder.h"
#include "scatterer/coreshell/coreshell.h"
#include "source/source.h"
#include "detector/detector.h"
#include "utils/numpy_interface.h"
#include "sets/sets.cpp"



class Experiment
{
    public:
        bool debug_mode = true;

        explicit Experiment(bool debug_mode = true) : debug_mode(debug_mode) {}


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
        void debug_print_state(const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const DetectorSet& detector_set) const {
            if (!debug_mode) return;
            debug_printf("----- Experiment Debug Info -----\n");
            debug_printf("SourceSet total combinations: %zu\n", source_set.total_combinations);
            debug_printf("ScattererSet total combinations: %zu\n", scatterer_set.total_combinations);
            debug_printf("DetectorSet total combinations: %zu\n", detector_set.total_combinations);
            debug_printf("SourceSet shape: ");
            for (size_t dim : source_set.shape) {
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

        template<double (BaseScatterer::*function)() const>
        pybind11::array_t<double> get_data(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const {

            if (debug_mode)
                this->debug_print_state(scatterer_set, source_set, detector_set);

            std::vector<size_t> array_shape;
            size_t total_iterations;

            if (detector_set.is_empty) {
                array_shape = concatenate_vector(source_set.shape, scatterer_set.shape);
                total_iterations = source_set.total_combinations * scatterer_set.total_combinations;
            }
            else {
                array_shape = concatenate_vector(source_set.shape, scatterer_set.shape, detector_set.shape);
                total_iterations = source_set.total_combinations * scatterer_set.total_combinations * detector_set.total_combinations;
            }

            debug_printf("get_scatterer_data: total_iterations = %zu\n", total_iterations);

            std::vector<double> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t flat_index = 0; flat_index < total_iterations; ++flat_index) {
                size_t idx; // Declare idx locally so each iteration has its own copy

                if (detector_set.is_empty) {
                    // 2D case: only source and scatterer
                    size_t i = flat_index / scatterer_set.total_combinations;
                    size_t j = flat_index % scatterer_set.total_combinations;
                    BaseSource source = source_set.get_source_by_index(i);

                    std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index(j, source);

                    idx = flatten_multi_index(array_shape, source.indices, scatterer_ptr->indices);

                    output_array[idx] = std::invoke(function, *scatterer_ptr);
                } else {
                    // 3D case: source, scatterer, and detector
                    size_t i = flat_index / (scatterer_set.total_combinations * detector_set.total_combinations);
                    size_t j = (flat_index / detector_set.total_combinations) % scatterer_set.total_combinations;
                    size_t k = flat_index % detector_set.total_combinations;
                    BaseSource source = source_set.get_source_by_index(i);

                    std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index(j, source);

                    Detector detector = detector_set.get_detector_by_index(k);
                    idx = flatten_multi_index(array_shape, source.indices, scatterer_ptr->indices, detector.indices);

                    output_array[idx] = std::invoke(function, *scatterer_ptr);
                }
            }

            debug_printf("get_scatterer_data: finished computation\n");
            return _vector_to_numpy(output_array, array_shape);
        }


        template<double (BaseScatterer::*function)() const > pybind11::array_t<double>
        get_data_sequential(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const {

            if (debug_mode)
                this->debug_print_state(scatterer_set, source_set, detector_set);

            std::vector<size_t> array_shape = {source_set.wavelength.size()};
            size_t full_size = source_set.wavelength.size();
            scatterer_set.validate_sequential_data(full_size);
            source_set.validate_sequential_data(full_size);
            debug_printf("get_scatterer_data_sequential: full_size = %zu\n", full_size);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {
                BaseSource source = source_set.get_source_by_index_sequential(idx);

                std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index_sequential(idx, source);

                output_array[idx] = std::invoke(function, *scatterer_ptr);

            }
            debug_printf("get_scatterer_data_sequential: finished computation\n");
            return _vector_to_numpy(output_array, {full_size});
        }

        pybind11::array_t<double> get_coupling(
            const ScattererSet& scatterer_set,
            const BaseSourceSet &source_set,
            const DetectorSet &detector_set) const {

            if (debug_mode)
                this->debug_print_state(scatterer_set, source_set, detector_set);

            std::vector<size_t> array_shape = concatenate_vector(source_set.shape, scatterer_set.shape, detector_set.shape);
            size_t total_iterations = source_set.total_combinations * scatterer_set.total_combinations * detector_set.total_combinations;
            debug_printf("get_scatterer_coupling: total_iterations = %zu\n", total_iterations);

            std::vector<double> output_array(total_iterations);

            #pragma omp parallel for
            for (size_t idx_flat = 0; idx_flat < total_iterations; ++idx_flat) {
                size_t i = idx_flat / (scatterer_set.total_combinations * detector_set.total_combinations);
                size_t j = (idx_flat / detector_set.total_combinations) % scatterer_set.total_combinations;
                size_t k = idx_flat % detector_set.total_combinations;

                BaseSource source = source_set.get_source_by_index(i);

                Detector detector = detector_set.get_detector_by_index(k);

                std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index(j, source);

                detector.medium_refractive_index = scatterer_ptr->medium_refractive_index;

                size_t idx = flatten_multi_index(array_shape, source.indices, scatterer_ptr->indices, detector.indices);
                output_array[idx] = detector.get_coupling(*scatterer_ptr);
            }
            debug_printf("get_scatterer_coupling: finished computation\n");
            return _vector_to_numpy(output_array, array_shape);
        }

        pybind11::array_t<double>
        get_coupling_sequential(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const {

            std::vector<size_t> array_shape = {source_set.wavelength.size()};
            size_t full_size = source_set.wavelength.size();
            scatterer_set.validate_sequential_data(full_size);
            source_set.validate_sequential_data(full_size);
            detector_set.validate_sequential_data(full_size);

            std::vector<double> output_array(full_size);

            #pragma omp parallel for
            for (size_t idx = 0; idx < full_size; ++idx) {

                BaseSource source = source_set.get_source_by_index_sequential(idx);

                std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index_sequential(idx, source);

                Detector detector = detector_set.get_detector_by_index_sequential(idx);

                detector.medium_refractive_index = scatterer_ptr->medium_refractive_index;

                output_array[idx] = detector.get_coupling(*scatterer_ptr);

            }

            return _vector_to_numpy(output_array, {full_size});
        }

};
