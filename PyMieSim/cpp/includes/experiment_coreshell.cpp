#pragma once

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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Ci][w],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n],
            source,
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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_material[Si][w],
            coreshellSet.n_medium[n],
            source,
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

    SOURCE::Planewave source = SOURCE::Planewave(
        sourceSet.wavelength[w],
        sourceSet.jones_vector[j],
        sourceSet.amplitude[w]
    );

    CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
        coreshellSet.core_diameter[Cd],
        coreshellSet.shell_width[Sd],
        coreshellSet.core_material[Ci][w],
        coreshellSet.shell_material[Si][w],
        coreshellSet.n_medium[n],
        source,
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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n],
            source,
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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Cm][w],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n],
            source
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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_material[Sm][w],
            coreshellSet.n_medium[n],
            source
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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Cm][w],
            coreshellSet.shell_material[Sm][w],
            coreshellSet.n_medium[n],
            source
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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n],
            source
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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

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

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n],
            source
        );

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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

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

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Ci][w],
            coreshellSet.shell_index[Si],
            coreshellSet.n_medium[n],
            source
        );

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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

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

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_index[Ci],
            coreshellSet.shell_material[Si][w],
            coreshellSet.n_medium[n],
            source
        );

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

        SOURCE::Planewave source = SOURCE::Planewave(
            sourceSet.wavelength[w],
            sourceSet.jones_vector[j],
            sourceSet.amplitude[w]
        );

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

        CORESHELL::Scatterer scatterer = CORESHELL::Scatterer(
            coreshellSet.core_diameter[Cd],
            coreshellSet.shell_width[Sd],
            coreshellSet.core_material[Ci][w],
            coreshellSet.shell_material[Si][w],
            coreshellSet.n_medium[n],
            source
        );

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}


