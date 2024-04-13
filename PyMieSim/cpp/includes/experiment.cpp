#pragma once

#include "experiment.h"


//--------------------------------------SPHERE------------------------------------
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        SPHERE::State scatterer_state = SPHERE::State(
            sphereSet.diameter[d],
            sphereSet.material[i][w],
            sphereSet.n_medium[n]
        );

        SPHERE::Scatterer scatterer = SPHERE::Scatterer(scatterer_state, source_state, max_order+1);

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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        SPHERE::State scatterer_state = SPHERE::State(
            sphereSet.diameter[d],
            sphereSet.index[i],
            sphereSet.n_medium[n]
        );


        SPHERE::Scatterer scatterer = SPHERE::Scatterer(
            scatterer_state,
            source_state,
            max_order+1
        );

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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        SPHERE::State scatterer_state = SPHERE::State(
            sphereSet.diameter[d],
            sphereSet.material[i][w],
            sphereSet.n_medium[n]
        );


        SPHERE::Scatterer scatterer = SPHERE::Scatterer(
            scatterer_state,
            source_state
        );

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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        SPHERE::State scatterer_state = SPHERE::State(
            sphereSet.diameter[d],
            sphereSet.index[i],
            sphereSet.n_medium[n]
        );


        SPHERE::Scatterer Scat = SPHERE::Scatterer(
            scatterer_state,
            source_state,
            max_order
        );

        output_array[idx] = (Scat.*function)();
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        SPHERE::State scatterer_state = SPHERE::State(
            sphereSet.diameter[d],
            sphereSet.material[i][w],
            sphereSet.n_medium[n]
        );

        DETECTOR::State detector_state  = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        SPHERE::Scatterer scatterer = SPHERE::Scatterer(scatterer_state, source_state);

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        SPHERE::State scatterer_state = SPHERE::State(
            sphereSet.diameter[d],
            sphereSet.index[i],
            sphereSet.n_medium[n]
        );


        DETECTOR::State detector_state = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        SPHERE::Scatterer scatterer = SPHERE::Scatterer(
            scatterer_state,
            source_state
        );

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

        double coupling = detector.get_coupling(scatterer);

        output_array[idx] = abs(coupling);
    }

    return vector_to_numpy(output_array, array_shape);
}


//--------------------------------------CYLINDER------------------------------------
template<typename Function>
pybind11::array_t<complex128> Experiment::get_cylinder_coefficient(Function function, size_t max_order) const
{
    if (cylinderSet.is_material)
        return get_cylinder_coefficient_material(function, max_order);
    else
        return get_cylinder_coefficient_index(function, max_order);
}

template<typename Function>
pybind11::array_t<double> Experiment::get_cylinder_data(Function function, size_t max_order) const
{
    if (cylinderSet.is_material)
        return get_cylinder_data_material(function, max_order);
    else
        return get_cylinder_data_index(function, max_order);
}

pybind11::array_t<double> Experiment::get_cylinder_coupling() const
{
    if (cylinderSet.is_material)
        return get_cylinder_coupling_bound();
    else
        return get_cylinder_coupling_unbound();
}


template<typename Function>
pybind11::array_t<complex128> Experiment::get_cylinder_coefficient_material(Function function, size_t max_order) const
{
    using namespace CYLINDER;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        cylinderSet.shape
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CYLINDER::State scatterer_state = CYLINDER::State(
            cylinderSet.diameter[d],
            cylinderSet.material[i][w],
            cylinderSet.n_medium[n]
        );

        CYLINDER::Scatterer Scat = CYLINDER::Scatterer(
            scatterer_state,
            source_state,
            max_order + 1
        );

        output_array[idx] = (Scat.*function)()[max_order];
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<complex128> Experiment::get_cylinder_coefficient_index(Function function, size_t max_order) const
{
    using namespace CYLINDER;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        cylinderSet.shape
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CYLINDER::State scatterer_state = CYLINDER::State(
            cylinderSet.diameter[d],
            cylinderSet.index[i],
            cylinderSet.n_medium[n]
        );

        CYLINDER::Scatterer Scat = CYLINDER::Scatterer(
            scatterer_state,
            source_state,
            max_order + 1
        );

        output_array[idx] = (Scat.*function)()[max_order];
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<double> Experiment::get_cylinder_data_material(Function function, size_t max_order) const
{
    using namespace CYLINDER;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        cylinderSet.shape
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CYLINDER::State scatterer_state = CYLINDER::State(
            cylinderSet.diameter[d],
            cylinderSet.material[i][w],
            cylinderSet.n_medium[n]
        );

        CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
            scatterer_state,
            source_state,
            max_order
        );

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<double> Experiment::get_cylinder_data_index(Function function, size_t max_order) const
{
    using namespace CYLINDER;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        cylinderSet.shape
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CYLINDER::State scatterer_state = CYLINDER::State(
            cylinderSet.diameter[d],
            cylinderSet.index[i],
            cylinderSet.n_medium[n]
        );

        CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
            scatterer_state,
            source_state,
            max_order
        );

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


pybind11::array_t<double> Experiment::get_cylinder_coupling_bound() const
{
    using namespace CYLINDER;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        cylinderSet.shape,
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CYLINDER::State scatterer_state = CYLINDER::State(
            cylinderSet.diameter[d],
            cylinderSet.material[i][w],
            cylinderSet.n_medium[n]
        );

        DETECTOR::State detector_state = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
            scatterer_state,
            source_state
        );

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}






pybind11::array_t<double> Experiment::get_cylinder_coupling_unbound() const
{
    using namespace CYLINDER;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        cylinderSet.shape,
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

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CYLINDER::State scatterer_state = CYLINDER::State(
            cylinderSet.diameter[d],
            cylinderSet.index[i],
            cylinderSet.n_medium[n]
        );

        DETECTOR::State detector_state = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        CYLINDER::Scatterer scatterer = CYLINDER::Scatterer(
            scatterer_state,
            source_state
        );

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

        output_array[idx] = abs(detector.get_coupling(scatterer));
    }

    return vector_to_numpy(output_array, array_shape);
}




//--------------------------------------CORESHELL------------------------------------
template<typename Function>
pybind11::array_t<complex128> Experiment::get_coreshell_coefficient(Function function, size_t max_order) const
{
    if (coreshellSet.core_is_material && coreshellSet.shell_is_material)
        return get_coreshell_coefficient_core_material_shell_material(function, max_order);

    if (coreshellSet.core_is_material && !coreshellSet.shell_is_material)
        return get_coreshell_coefficient_core_material_shell_index(function, max_order);

    if (!coreshellSet.core_is_material && coreshellSet.shell_is_material)
        return get_coreshell_coefficient_core_index_shell_material(function, max_order);

    if (!coreshellSet.core_is_material && !coreshellSet.shell_is_material)
        return get_coreshell_coefficient_core_index_shell_index(function, max_order);
}

template<typename Function>
pybind11::array_t<double> Experiment::get_coreshell_data(Function function, size_t max_order) const
{
    if (coreshellSet.core_is_material && coreshellSet.shell_is_material)
        return get_coreshell_data_core_material_shell_material(function, max_order);

    if (coreshellSet.core_is_material && !coreshellSet.shell_is_material)
        return get_coreshell_data_core_material_shell_index(function, max_order);

    if (!coreshellSet.core_is_material && coreshellSet.shell_is_material)
        return get_coreshell_data_core_index_shell_material(function, max_order);

    if (!coreshellSet.core_is_material && !coreshellSet.shell_is_material)
        return get_coreshell_data_core_index_shell_index(function, max_order);
}


pybind11::array_t<double> Experiment::get_coreshell_coupling() const
{
    if (coreshellSet.core_is_material && coreshellSet.shell_is_material)
        return get_coreshell_coupling_core_material_shell_material();

    if (coreshellSet.core_is_material && !coreshellSet.shell_is_material)
        return get_coreshell_coupling_core_material_shell_index();

    if (!coreshellSet.core_is_material && coreshellSet.shell_is_material)
        return get_coreshell_coupling_core_index_shell_material();

    if (!coreshellSet.core_is_material && !coreshellSet.shell_is_material)
        return get_coreshell_coupling_core_index_shell_index();
}


template<typename Function>
pybind11::array_t<complex128> Experiment::get_coreshell_coefficient_core_material_shell_index(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<complex128> output_array(full_size);

    #pragma omp parallel for collapse(7)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Ci][w],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state,
            max_order + 1
        );

        output_array[idx] = (scatterer.*function)()[max_order];
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<complex128> Experiment::get_coreshell_coefficient_core_index_shell_material(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<complex128> output_array(full_size);

    #pragma omp parallel for collapse(7)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_material[Si][w],
            coreshellSet.n_medium[n]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state,
            max_order + 1
        );

        output_array[idx] = (scatterer.*function)()[max_order];
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<complex128> Experiment::get_coreshell_coefficient_core_material_shell_material(Function function, size_t max_order) const
{
using namespace CORESHELL;

std::vector<size_t> array_shape = concatenate_vector(
    sourceSet.shape,
    coreshellSet.shape
);

size_t full_size = get_vector_sigma(array_shape);

std::vector<complex128> output_array(full_size);

#pragma omp parallel for collapse(7)
for (size_t w=0; w<array_shape[0]; ++w)
for (size_t j=0; j<array_shape[1]; ++j)
for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
for (size_t Si=0; Si<array_shape[5]; ++Si)
for (size_t n=0; n<array_shape[6]; ++n)
{
    size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

    SOURCE::State source_state = SOURCE::State(
        sourceSet.wavelength[w],
        sourceSet.jones_vector[j],
        sourceSet.amplitude[w]
    );

    CORESHELL::State scatterer_state = CORESHELL::State(
        coreshellSet.core_diameter[Cd],
        coreshellSet.shell_width[Sd],
        coreshellSet.core_material[Ci][w],
        coreshellSet.shell_material[Si][w],
        coreshellSet.n_medium[n]
    );

    CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
        scatterer_state,
        source_state,
        max_order + 1
    );

    output_array[idx] = (scatterer.*function)()[max_order];
}

return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<complex128> Experiment::get_coreshell_coefficient_core_index_shell_index(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<complex128> output_array(full_size);

    #pragma omp parallel for collapse(7)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state,
            max_order + 1
        );

        output_array[idx] = (scatterer.*function)()[max_order];
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<double> Experiment::get_coreshell_data_core_material_shell_index(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(7)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Cm=0; Cm<array_shape[4]; ++Cm)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Cm, Si, n}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Cm][w],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state,
            max_order
        );

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<double> Experiment::get_coreshell_data_core_index_shell_material(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(7)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Sm=0; Sm<array_shape[5]; ++Sm)
    for (size_t n=0; n<array_shape[6]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Sm, n}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_material[Sm][w],
            coreshellSet.n_medium[n]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state,
            max_order
        );

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<double> Experiment::get_coreshell_data_core_material_shell_material(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(7)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Cm=0; Cm<array_shape[4]; ++Cm)
    for (size_t Sm=0; Sm<array_shape[5]; ++Sm)
    for (size_t n=0; n<array_shape[6]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Cm, Sm, n}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Cm][w],
            coreshellSet.shell_material[Sm][w],
            coreshellSet.n_medium[n]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state,
            max_order
        );

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


template<typename Function>
pybind11::array_t<double> Experiment::get_coreshell_data_core_index_shell_index(Function function, size_t max_order) const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(7)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state,
            max_order
        );

        output_array[idx] = (scatterer.*function)();
    }

    return vector_to_numpy(output_array, array_shape);
}


pybind11::array_t<double> Experiment::get_coreshell_coupling_core_index_shell_index() const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);


    #pragma omp parallel for collapse(12)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n]
        );

        DETECTOR::State detector_state = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state
        );

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}




pybind11::array_t<double> Experiment::get_coreshell_coupling_core_material_shell_index() const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(12)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Ci][w],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n]
        );

        DETECTOR::State detector_state = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state
        );

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}




pybind11::array_t<double> Experiment::get_coreshell_coupling_core_index_shell_material() const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(12)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_material[Si][w],
            coreshellSet.n_medium[n]
        );

        DETECTOR::State detector_state = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state
        );

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}


pybind11::array_t<double> Experiment::get_coreshell_coupling_core_material_shell_material() const
{
    using namespace CORESHELL;

    std::vector<size_t> array_shape = concatenate_vector(
        sourceSet.shape,
        coreshellSet.shape,
        detectorSet.shape
    );

    size_t full_size = get_vector_sigma(array_shape);

    std::vector<double> output_array(full_size);

    #pragma omp parallel for collapse(12)
    for (size_t w=0; w<array_shape[0]; ++w)
    for (size_t j=0; j<array_shape[1]; ++j)
    for (size_t Cd=0; Cd<array_shape[2]; ++Cd)
    for (size_t Sd=0; Sd<array_shape[3]; ++Sd)
    for (size_t Ci=0; Ci<array_shape[4]; ++Ci)
    for (size_t Si=0; Si<array_shape[5]; ++Si)
    for (size_t n=0; n<array_shape[6]; ++n)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({w, j, Cd, Sd, Ci, Si, n, s, na, p, g, f}, array_shape);

        SOURCE::State source_state = SOURCE::State(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::State scatterer_state = CORESHELL::State(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Ci][w],
            coreshellSet.shell_material[Si][w],
            coreshellSet.n_medium[n]
        );

        DETECTOR::State  detector_state = DETECTOR::State(
            detectorSet.scalar_fields[s],
            detectorSet.NA[na],
            detectorSet.phi_offset[p],
            detectorSet.gamma_offset[g],
            detectorSet.polarization_filter[f],
            detectorSet.rotation_angle[s],
            detectorSet.coherent,
            detectorSet.point_coupling
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            scatterer_state,
            source_state
        );

        DETECTOR::Detector detector = DETECTOR::Detector(detector_state);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}


