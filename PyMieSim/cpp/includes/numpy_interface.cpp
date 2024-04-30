#pragma once

#include "definitions.cpp"

template <typename T>
std::vector<size_t>
get_stride(std::vector<size_t> dimension)
{
    std::reverse(dimension.begin(), dimension.end());

    std::vector<size_t> stride;
    stride.push_back( sizeof(T) );

    for (size_t i=0; i<dimension.size()-1; ++i)
        stride.push_back( stride[i] * dimension[i] );

    std::reverse(stride.begin(), stride.end());

    return stride;

}

template <typename T>
pybind11::array_t<T>
vector_to_numpy(const std::vector<T> &vector, const std::vector<size_t> &dimension)
{
    pybind11::array_t<T>  numpy_array;

    std::vector<T> * output_vector = new std::vector<T>;
    (*output_vector) = vector;

    std::vector<size_t> stride = get_stride<T>(dimension);

    pybind11::capsule free_when_done(
        output_vector->data(), [](void *f) {T *foo = reinterpret_cast<T *>(f); delete []foo; }
    );

    numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        output_vector->data(),
        free_when_done
      );

    return numpy_array;
}


template <typename T>
pybind11::array_t<T>
vector_to_numpy(const std::vector<T> &vector)
{
    pybind11::array_t<T>  numpy_array;

    std::vector<T> * output_vector = new std::vector<T>;
    (*output_vector) = vector;

    std::vector<size_t> stride = {sizeof(T)};
    std::vector<size_t> dimension ={vector.size()};

    pybind11::capsule free_when_done(
        output_vector->data(), [](void *f) {T *foo = reinterpret_cast<T *>(f); delete []foo; }
    );

    numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        output_vector->data(),
        free_when_done
      );

    return numpy_array;
}


template<typename T>
pybind11::array_t<T>
vector_to_numpy(std::vector<T> &vector, const std::vector<size_t> &dimension, const std::vector<size_t> &stride){

    pybind11::capsule free_when_done(
        vector.data(), [](void *f) { T *foo = reinterpret_cast<T *>(f); }
    );

    pybind11::array_t<T> numpy_array = pybind11::array_t<T>(
        dimension,
        stride,
        vector.data(),
        free_when_done
      );

    return numpy_array;
}

template<typename T>
inline pybind11::array_t<T>
vector_to_numpy(std::vector<T>&& passthrough)
{
    auto* ptr = new std::vector<T>(std::move(passthrough));

    const pybind11::capsule freeWhenDone(
        ptr, [](void *toFree) { delete static_cast<std::vector<T> *>(toFree); }
    );

    auto numpy_array = pybind11::array_t<T>(
        {ptr->size()},
        {sizeof(T)},
        ptr->data(),
        freeWhenDone
    );

    return numpy_array;
}

template<typename T>
inline pybind11::array_t<T>
vector_to_numpy_copy(std::vector<T> passthrough)
{
    auto* ptr = new std::vector<T>(std::move(passthrough));

    const pybind11::capsule freeWhenDone(
        ptr, [](void *toFree) { delete static_cast<std::vector<T> *>(toFree); }
    );

    auto numpy_array = pybind11::array_t<T>(
        {ptr->size()},
        {sizeof(T)},
        ptr->data(),
        freeWhenDone
    );

    return numpy_array;
}


template<typename T>
inline pybind11::array_t<T>
vector_to_numpy_copy(std::vector<T> passthrough, const std::vector<size_t> &shape)
{
    auto* ptr = new std::vector<T>(std::move(passthrough));

    const pybind11::capsule freeWhenDone(
        ptr, [](void *toFree) { delete static_cast<std::vector<T> *>(toFree); }
    );

    auto numpy_array = pybind11::array_t<T>(
        {ptr->size()},
        {sizeof(T)},
        ptr->data(),
        freeWhenDone
    );

    numpy_array.resize(shape);

    return numpy_array;
}


template<typename T>
inline pybind11::array_t<T>
vector_to_numpy(std::vector<T>&& passthrough, const std::vector<size_t> &shape)
{
    auto* ptr = new std::vector<T>(std::move(passthrough));

    const pybind11::capsule freeWhenDone(
        ptr, [](void *toFree) { delete static_cast<std::vector<T> *>(toFree); }
    );

    auto numpy_array = pybind11::array_t<T>(
        {ptr->size()},
        {sizeof(T)},
        ptr->data(),
        freeWhenDone
    );

    numpy_array.resize(shape);

    return numpy_array;
}

