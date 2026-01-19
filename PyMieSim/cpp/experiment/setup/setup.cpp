#include "setup.h"



// Example: You can also add a class-wide debug_print method.
void Setup::debug_print_state(const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const BaseDetectorSet& detector_set) const {
    if (!debug_mode) return;
    debug_printf("----- Setup Debug Info -----\n");
    debug_printf("SourceSet total combinations: %zu\n", source_set.total_combinations);
    debug_printf("ScattererSet total combinations: %zu\n", scatterer_set.total_combinations);
    debug_printf("DetectorSet total combinations: %zu\n", detector_set.total_combinations);
    debug_printf("SourceSet shape: ");
    for (size_t dim : source_set.shape) {
        debug_printf("%zu ", dim);
    }
    debug_printf("\n");
    debug_printf("ScattererSet shape: ");
    for (size_t dim : scatterer_set.shape) {
        debug_printf("%zu ", dim);
    }
    debug_printf("\n");
    debug_printf("DetectorSet shape: ");
    for (size_t dim : detector_set.shape) {
        debug_printf("%zu ", dim);
    }
    debug_printf("\n---------------------------------\n");
}


std::tuple<std::vector<double>, std::vector<size_t>>
Setup::get_coupling(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const BaseDetectorSet &detector_set) const {

    if (debug_mode)
        this->debug_print_state(scatterer_set, source_set, detector_set);

    std::vector<size_t> array_shape = concatenate_vector(source_set.shape, scatterer_set.shape, detector_set.shape);
    size_t total_iterations = source_set.total_combinations * scatterer_set.total_combinations * detector_set.total_combinations;
    debug_printf("get_scatterer_coupling: total_iterations = %zu\n", total_iterations);

    std::vector<double> output_array(total_iterations);

    #pragma omp parallel for
    for (long long idx_flat = 0; idx_flat < static_cast<long long>(total_iterations); ++idx_flat) {
        size_t i = idx_flat / (scatterer_set.total_combinations * detector_set.total_combinations);
        size_t j = (idx_flat / detector_set.total_combinations) % scatterer_set.total_combinations;
        size_t k = idx_flat % detector_set.total_combinations;

        std::shared_ptr<BaseSource> source = source_set.get_source_by_index(i);

        std::shared_ptr<BaseDetector> detector = detector_set.get_detector_by_index(k);

        std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index(j, source);

        detector->scatterer_medium_refractive_index = scatterer_ptr->medium_refractive_index;
        size_t idx = flatten_multi_index(array_shape, source->indices, scatterer_ptr->indices, detector->indices);
        output_array[idx] = detector->get_coupling(*scatterer_ptr);
    }
    debug_printf("get_scatterer_coupling: finished computation\n");

    return std::make_tuple(std::move(output_array), std::move(array_shape));
}

std::vector<double>
Setup::get_coupling_sequential(const ScattererSet& scatterer_set, const BaseSourceSet &source_set, const BaseDetectorSet &detector_set) const {

    std::vector<size_t> array_shape = {source_set.wavelength.size()};
    size_t full_size = source_set.wavelength.size();
    scatterer_set.validate_sequential_data(full_size);
    source_set.validate_sequential_data(full_size);
    detector_set.validate_sequential_data(full_size);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for
    for (long long idx = 0; idx < static_cast<long long>(full_size); ++idx) {

        std::shared_ptr<BaseSource> source = source_set.get_source_by_index_sequential(idx);

        std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index_sequential(idx, source);

        std::shared_ptr<BaseDetector> detector = detector_set.get_detector_by_index_sequential(idx);

        detector->scatterer_medium_refractive_index = scatterer_ptr->medium_refractive_index;

        output_array[idx] = detector->get_coupling(*scatterer_ptr);
    }

    return output_array;
}


std::tuple<std::vector<complex128>, std::vector<size_t>>
Setup::get_farfields(const ScattererSet& scatterer_set, const BaseSourceSet& source_set, const FibonacciMesh& mesh, const double distance) const
{
    // Head shape for indexing the source and scatterer grid only
    std::vector<size_t> head_shape = concatenate_vector(source_set.shape, scatterer_set.shape);

    // Full output shape: [..., 2, sampling]
    std::vector<size_t> array_shape = head_shape;
    array_shape.push_back(2);                       // channel: 0 -> theta, 1 -> phi
    array_shape.push_back(mesh.sampling);       // samples along the mesh

    const size_t total_iterations = source_set.total_combinations * scatterer_set.total_combinations;
    const size_t stride_channel   = mesh.sampling;      // size of one trace
    const size_t stride_block     = 2 * stride_channel;     // theta + phi per combo

    // Allocate the full array
    size_t total_size = 1;
    for (size_t d : array_shape)
        total_size *= d;

    std::vector<complex128> output_array(total_size, 0.0);

    debug_printf("get_scatterer_farfields: total_iterations = %zu\n", total_iterations);

    #pragma omp parallel for
    for (long long idx_flat = 0; idx_flat < static_cast<long long>(total_iterations); ++idx_flat) {
        size_t i = idx_flat / scatterer_set.total_combinations;  // source index
        size_t j = idx_flat % scatterer_set.total_combinations;  // scatterer index

        std::shared_ptr<BaseSource> source = source_set.get_source_by_index(i);
        std::unique_ptr<BaseScatterer> scatterer_ptr = scatterer_set.get_scatterer_ptr_by_index(j, std::move(source));

        // Linear index over the head shape only
        size_t idx_head = flatten_multi_index(head_shape, source->indices, scatterer_ptr->indices);

        // Compute fields
        auto [phi_field, theta_field] = scatterer_ptr->compute_unstructured_farfields(mesh, distance);

        // Sanity check in debug builds
        if (phi_field.size() != stride_channel || theta_field.size() != stride_channel) {
            // handle error as you prefer
            continue;
        }

        // Write into the last two axes
        size_t base = idx_head * stride_block;

        std::copy(theta_field.begin(), theta_field.end(), output_array.begin() + base + 0 * stride_channel);

        std::copy(phi_field.begin(),   phi_field.end(), output_array.begin() + base + 1 * stride_channel);
    }

    debug_printf("get_scatterer_farfields: finished computation\n");
    return std::make_tuple(std::move(output_array), std::move(array_shape));
}
