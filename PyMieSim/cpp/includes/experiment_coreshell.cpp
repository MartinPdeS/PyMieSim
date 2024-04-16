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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
for (size_t wl=0; wl<array_shape[0]; ++wl)
for (size_t jv=0; jv<array_shape[1]; ++jv)
for (size_t cd=0; cd<array_shape[2]; ++cd)
for (size_t sw=0; sw<array_shape[3]; ++sw)
for (size_t cm=0; cm<array_shape[4]; ++cm)
for (size_t sm=0; sm<array_shape[5]; ++sm)
for (size_t mi=0; mi<array_shape[6]; ++mi)
{
    size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

    SOURCE::Planewave source = sourceSet.to_object(wl, jv);

    CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi, s, na, p, g, f}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

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

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi, s, na, p, g, f}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

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

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi, s, na, p, g, f}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

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

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

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
    for (size_t wl=0; wl<array_shape[0]; ++wl)
    for (size_t jv=0; jv<array_shape[1]; ++jv)
    for (size_t cd=0; cd<array_shape[2]; ++cd)
    for (size_t sw=0; sw<array_shape[3]; ++sw)
    for (size_t cm=0; cm<array_shape[4]; ++cm)
    for (size_t sm=0; sm<array_shape[5]; ++sm)
    for (size_t mi=0; mi<array_shape[6]; ++mi)
    for (size_t s=0; s<array_shape[7]; ++s)
    for (size_t na=0; na<array_shape[8]; ++na)
    for (size_t p=0; p<array_shape[9]; ++p)
    for (size_t g=0; g<array_shape[10]; ++g)
    for (size_t f=0; f<array_shape[11]; ++f)
    {
        size_t idx = flatten_multi_index({wl, jv, cd, sw, cm, sm, mi, s, na, p, g, f}, array_shape);

        SOURCE::Planewave source = sourceSet.to_object(wl, jv);

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

        CORESHELL::Scatterer scatterer = coreshellSet.to_object(wl, cd, sw, cm, sm, mi, source);

        output_array[idx] = abs( detector.get_coupling(scatterer) );
    }

    return vector_to_numpy(output_array, array_shape);
}


