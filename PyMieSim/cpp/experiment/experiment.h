#pragma once

#include <cstdarg>    // For va_list, va_start, va_end
#include <cstdio>     // For printf, vprintf
#include "scatterer/sphere/sphere.h"
#include "scatterer/cylinder/cylinder.h"
#include "scatterer/coreshell/coreshell.h"
#include "source/source.h"
#include "detector/detector.h"
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

        /**
         * @brief Prints the current state of the experiment for debugging purposes.
         * @param scatterer_set The set of scatterers involved in the experiment.
         * @param source_set The set of sources involved in the experiment.
         * @param detector_set The set of detectors involved in the experiment.
        */
        void debug_print_state(const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const DetectorSet& detector_set) const;

        /**
         * @brief Flattens a multi-dimensional index into a single index.
         * @tparam MultiIndexVectors Variadic template for multiple index vectors.
         * @param dimensions The dimensions of the multi-dimensional array.
         * @param multi_indices The multi-dimensional indices to be flattened.
         * @return The flattened single index.
         */
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

        /**
         * @brief Concatenates multiple vectors into a single vector.
         * @tparam T The type of the vectors.
         * @tparam Ts Variadic template for multiple vector types.
         * @param first_vector The first vector to concatenate.
         * @param other_vectors The other vectors to concatenate.
         * @return The concatenated vector.
         */
        template <typename T, typename... Ts>
        T concatenate_vector(const T& first_vector, const Ts&... other_vectors) const
        {
            T output_vector = first_vector;
            (output_vector.insert(output_vector.end(), other_vectors.begin(), other_vectors.end()), ...);
            return output_vector;
        }

        /**
         * @brief Generic method to get scatterer data by invoking a member function.
         * @tparam function Pointer to the member function of BaseScatterer to invoke.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A tuple containing a numpy array of the requested data and the shape of the array.
         */
        template<double (BaseScatterer::*function)() const>
        std::tuple<std::vector<double>, std::vector<size_t>> get_data(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const {

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
            return std::make_tuple(std::move(output_array), std::move(array_shape));
        }


        /**
         * @brief Generic method to get scatterer data sequentially by invoking a member function.
         * @tparam function Pointer to the member function of BaseScatterer to invoke.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A numpy array containing the requested data.
         */
        template<double (BaseScatterer::*function)() const >
        std::vector<double>
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
            return output_array;
        }

        /**
         * @brief Computes the coupling coefficient for given scatterers, sources, and detectors.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A tuple containing a numpy array of coupling coefficients and the shape of the array.
         */
        std::tuple<std::vector<double>, std::vector<size_t>>
        get_coupling(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const;

        /**
         * @brief Computes the coupling coefficient sequentially for given scatterers, sources, and detectors.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param detector_set The set of detectors.
         * @return A numpy array of coupling coefficients.
         */
        std::vector<double>
        get_coupling_sequential(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const DetectorSet &detector_set) const;

        /**
         * @brief Computes the far-field patterns for given scatterers, sources, and a Fibonacci mesh.
         * @param scatterer_set The set of scatterers.
         * @param source_set The set of sources.
         * @param mesh The Fibonacci mesh for far-field sampling.
         * @param distance The distance at which to compute the far-fields (default is 1).
         * @return A tuple containing a numpy array of far-field patterns and the shape of the array.
         */
        std::tuple<std::vector<complex128>, std::vector<size_t>>
        get_farfields(const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const FibonacciMesh& mesh, const double distance = 1) const;

};
