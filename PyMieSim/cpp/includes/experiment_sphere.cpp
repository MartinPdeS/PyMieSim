#pragma once


template<typename Function>
pybind11::array_t<complex128> Experiment::get_sphere_coefficient(Function function, size_t max_order) const
{
    if (sphereSet.is_material)
        return get_sphere_coefficient_material(function, max_order);
    else
        return get_sphere_coefficient_index(function, max_order);
}

template<typename Function>
pybind11::array_t<double> Experiment::get_sphere_data(Function function, size_t max_order) const
{
    if (sphereSet.is_material)
        return get_sphere_data_material(function, max_order);
    else
        return get_sphere_data_index(function, max_order);
}


pybind11::array_t<double> Experiment::get_sphere_coupling() const
{
    if (sphereSet.is_material)
        return get_sphere_coupling_material();
    else
        return get_sphere_coupling_index();
}


size_t Experiment::flatten_multi_index(const std::vector<size_t> &multi_index, const std::vector<size_t> &dimension) const
{
    size_t
        flatten_index = 0,
        multiplier = 1,
        dimension_size,
        index;

    for (size_t index_number=0; index_number<multi_index.size(); index_number++)
    {
        index = multi_index[index_number];
        multiplier = 1;
        for (size_t dim_number=index_number+1; dim_number<dimension.size(); dim_number++)
        {
            dimension_size = dimension[dim_number];
            multiplier *= dimension_size;
        }

        flatten_index += index * multiplier;
    }

    return flatten_index;
}

template<typename Function>
pybind11::array_t<complex128> Experiment::get_sphere_coefficient_material(Function function, size_t max_order) const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<complex128> output_array(full_size);

    #pragma omp parallel for collapse(5)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t d=0; d<array_shape[2]; ++d)
    for (size_t i=0; i<array_shape[3]; ++i)
    for (size_t n=0; n<array_shape[4]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(w, j);

        SPHERE::Scatterer scatterer = sphereSet.to_object(d, i, w, n, source);

        output_array[idx] = (scatterer.*function)()[max_order];
    }

  return vector_to_numpy(output_array, array_shape);
}

template<typename Function>
pybind11::array_t<complex128> Experiment::get_sphere_coefficient_index(Function function, size_t max_order) const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<complex128> output_array(full_size);

    #pragma omp parallel for collapse(5)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t d=0; d<array_shape[2]; ++d)
    for (size_t i=0; i<array_shape[3]; ++i)
    for (size_t n=0; n<array_shape[4]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(w, j);

        SPHERE::Scatterer scatterer = sphereSet.to_object(d, i, w, n, source);

        output_array[idx] = (scatterer.*function)()[max_order];
    }

  return vector_to_numpy(output_array, array_shape);
}

template<typename Function>
pybind11::array_t<double> Experiment::get_sphere_data_material(Function function, size_t max_order) const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(5)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t d=0; d<array_shape[2]; ++d)
    for (size_t i=0; i<array_shape[3]; ++i)
    for (size_t n=0; n<array_shape[4]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(w, j);

        SPHERE::Scatterer scatterer = sphereSet.to_object(d, i, w, n, source);

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<double> Experiment::get_sphere_data_index(Function function, size_t max_order) const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(5)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t d=0; d<array_shape[2]; ++d)
    for (size_t i=0; i<array_shape[3]; ++i)
    for (size_t n=0; n<array_shape[4]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, d, i, n}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(w, j);

        SPHERE::Scatterer scatterer = sphereSet.to_object(d, i, w, n, source);

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}

pybind11::array_t<double> Experiment::get_sphere_coupling_material() const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(10)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t d=0; d<array_shape[2]; ++d)
    for (size_t i=0; i<array_shape[3]; ++i)
    for (size_t n=0; n<array_shape[4]; ++n)
    for (size_t s=0; s<array_shape[5]; ++s)
    for (size_t na=0; na<array_shape[6]; ++na)
    for (size_t p=0; p<array_shape[7]; ++p)
    for (size_t g=0; g<array_shape[8]; ++g)
    for (size_t f=0; f<array_shape[9]; ++f)
    {
        size_t idx = flatten_multi_index({w, j, d, i, n, s, na, p, g, f}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(w, j);

        DETECTOR::Detector detector = DETECTOR::Detector(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        SPHERE::Scatterer scatterer = sphereSet.to_object(d, i, w, n, source);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}



pybind11::array_t<double> Experiment::get_sphere_coupling_index() const
{
    using namespace SPHERE;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        sphereSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(10)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t d=0; d<array_shape[2]; ++d)
    for (size_t i=0; i<array_shape[3]; ++i)
    for (size_t n=0; n<array_shape[4]; ++n)
    for (size_t s=0; s<array_shape[5]; ++s)
    for (size_t na=0; na<array_shape[6]; ++na)
    for (size_t p=0; p<array_shape[7]; ++p)
    for (size_t g=0; g<array_shape[8]; ++g)
    for (size_t f=0; f<array_shape[9]; ++f)
    {
        size_t idx = flatten_multi_index({w, j, d, i, n, s, na, p, g, f}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(w, j);

        DETECTOR::Detector detector = DETECTOR::Detector(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        SPHERE::Scatterer scatterer = sphereSet.to_object(d, i, w, n, source);

        double coupling = detector.get_coupling(scatterer);

        output_array[idx] = abs(coupling);
    }

    return vector_to_numpy(output_array, array_shape);
}
