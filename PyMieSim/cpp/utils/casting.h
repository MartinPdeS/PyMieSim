#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <sstream>
#include <complex>
#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/complex.h>
#include <pybind11/stl.h>
#include <experiment/polarization_set/polarization_set.h>

namespace py = pybind11;
using complex128 = std::complex<double>;

namespace Casting {


    template <typename dtype>
    std::vector<dtype> cast_py_to_vector(
        const py::object& object,
        const std::string& units = ""
    ) {
        if (!units.empty()) {
            try {
                py::object magnitude = py::reinterpret_borrow<py::object>(
                    object.attr("to")(units).attr("magnitude")
                );

                return cast_py_to_vector<dtype>(magnitude);
            }
            catch (const py::error_already_set&) {
                std::ostringstream oss;
                oss << "Failed to convert object to '" << units
                    << "'. Ensure the object has compatible units using PyMieSim.units.ureg.";
                throw std::invalid_argument(oss.str());
            }
        }

        if (py::isinstance<py::array>(object)) {
            py::array array = py::reinterpret_borrow<py::array>(object);

            if (array.ndim() == 0) {
                return {py::cast<dtype>(array.attr("item")())};
            }

            if (array.ndim() == 1) {
                return object.cast<std::vector<dtype>>();
            }

            throw std::invalid_argument("Only scalar or 1D array inputs are supported.");
        }

        try {
            return {object.cast<dtype>()};
        }
        catch (const py::cast_error&) {
        }

        if (py::isinstance<py::sequence>(object) && !py::isinstance<py::str>(object)) {
            return object.cast<std::vector<dtype>>();
        }

        throw std::invalid_argument("Object cannot be converted to a scalar or sequence.");
    }

    template <typename dtype>
    std::vector<dtype> cast_py_to_broadcasted_vector(
        const std::string& name,
        const py::object& object,
        const size_t target_size,
        const std::string& units = ""
    ) {
        if (target_size == 0) {
            std::ostringstream oss;
            oss << "Parameter '" << name << "' cannot be broadcast because target_size is 0.";
            throw std::invalid_argument(oss.str());
        }

        std::vector<dtype> values = cast_py_to_vector<dtype>(object, units);

        if (values.empty()) {
            std::ostringstream oss;
            oss << "Parameter '" << name << "' is empty. Provide a scalar or a non empty array.";
            throw std::invalid_argument(oss.str());
        }

        if (values.size() == 1) {
            return std::vector<dtype>(target_size, values[0]);
        }

        if (values.size() != target_size) {
            std::ostringstream oss;
            oss << "Inconsistent sizes: '" << name << "' has size " << values.size() << " but expected 1 or " << target_size << ".";
            throw std::invalid_argument(oss.str());
        }

        return values;
    }

    template <typename MaterialSetType, typename RefractiveIndexType, typename BaseClass, typename ConstantClass>
    MaterialSetType create_material_set_from_pyobject(
        const py::object& material_object,
        const std::string& material_name
    ) {
        if (py::isinstance<MaterialSetType>(material_object)) {
            return py::cast<MaterialSetType>(material_object);
        }

        if (py::isinstance<BaseClass>(material_object)) {
            return MaterialSetType(
                std::vector<std::shared_ptr<BaseClass>>{
                    py::cast<std::shared_ptr<BaseClass>>(material_object)
                }
            );
        }

        if (
            py::isinstance<py::float_>(material_object) ||
            py::isinstance<py::int_>(material_object) ||
            PyComplex_Check(material_object.ptr())
        ) {
            return MaterialSetType(
                std::vector<std::shared_ptr<BaseClass>>{
                    std::make_shared<ConstantClass>(
                        material_object.cast<RefractiveIndexType>()
                    )
                }
            );
        }

        if (py::isinstance<py::sequence>(material_object) && !py::isinstance<py::str>(material_object)) {
            py::sequence material_sequence = py::reinterpret_borrow<py::sequence>(material_object);

            std::vector<std::shared_ptr<BaseClass>> materials;
            materials.reserve(py::len(material_sequence));

            for (size_t i = 0; i < py::len(material_sequence); ++i) {
                py::object item = py::reinterpret_borrow<py::object>(material_sequence[i]);

                if (py::isinstance<BaseClass>(item)) {
                    materials.push_back(py::cast<std::shared_ptr<BaseClass>>(item));
                    continue;
                }

                if (
                    py::isinstance<py::float_>(item) ||
                    py::isinstance<py::int_>(item) ||
                    PyComplex_Check(item.ptr())
                ) {
                    materials.push_back(
                        std::make_shared<ConstantClass>(
                            item.cast<RefractiveIndexType>()
                        )
                    );
                    continue;
                }

                throw std::runtime_error(
                    "Invalid item in " + material_name + " sequence at index " + std::to_string(i) +
                    ". Expected a scalar refractive index or a Material instance."
                );
            }

            return MaterialSetType(materials);
        }

        throw std::runtime_error(
            "Invalid type for " + material_name + ". Expected a MaterialSet, "
            "a scalar refractive index, a sequence of refractive indices, "
            "a Material instance, or a sequence mixing refractive indices and Material instances."
        );
    }


    inline PolarizationSet cast_py_to_polarization_set(
        const py::object& polarization,
        const size_t target_size
    )
    {
        if (py::isinstance<PolarizationSet>(polarization))
            return polarization.cast<PolarizationSet>();

        if (polarization.is_none())
            return PolarizationSet(PolarizationState());

        if (py::isinstance<PolarizationState>(polarization))
            return PolarizationSet(
                std::vector<PolarizationState>(
                    target_size,
                    polarization.cast<PolarizationState>()
                )
            );


        try {
            return PolarizationSet(
                Casting::cast_py_to_broadcasted_vector<double>(
                    "polarization",
                    polarization,
                    target_size,
                    "radian"
                )
            );
        }
        catch (const std::exception&) {
        }

        if (py::isinstance<py::list>(polarization) || py::isinstance<py::tuple>(polarization)) {
            py::sequence sequence = py::reinterpret_borrow<py::sequence>(polarization);

            if (py::len(sequence) == 0) {
                throw std::invalid_argument(
                    "Invalid type for polarization. Empty sequences are not allowed."
                );
            }

            bool all_are_polarization_state = true;

            for (size_t index = 0; index < py::len(sequence); ++index) {
                py::object item = py::reinterpret_borrow<py::object>(sequence[index]);

                if (!py::isinstance<PolarizationState>(item)) {
                    all_are_polarization_state = false;
                    break;
                }
            }

            if (all_are_polarization_state) {
                std::vector<PolarizationState> state_vector;
                state_vector.reserve(py::len(sequence));

                for (size_t index = 0; index < py::len(sequence); ++index) {
                    py::object item = py::reinterpret_borrow<py::object>(sequence[index]);
                    state_vector.push_back(item.cast<PolarizationState>());
                }

                if (state_vector.size() == 1) {
                    state_vector = std::vector<PolarizationState>(target_size, state_vector[0]);
                }
                else if (state_vector.size() != target_size) {
                    throw std::invalid_argument(
                        "Inconsistent sizes: 'polarization' has size " +
                        std::to_string(state_vector.size()) +
                        " but expected 1 or " +
                        std::to_string(target_size) + "."
                    );
                }

                return PolarizationSet(state_vector);
            }
        }

        throw std::invalid_argument(
            "Invalid type for polarization. Expected None, a scalar angle, a sequence of angles, "
            "a PolarizationState, or a sequence of PolarizationState."
        );
    }

    inline PolarizationSet cast_py_to_polarization_set(
        const py::object& polarization
    )
    {

        if (py::isinstance<PolarizationSet>(polarization))
            return polarization.cast<PolarizationSet>();


        if (polarization.is_none())
            return PolarizationSet(PolarizationState());


        if (py::isinstance<PolarizationState>(polarization))
            return PolarizationSet(
                std::vector<PolarizationState>{
                    polarization.cast<PolarizationState>()
                }
            );


        try {
            return PolarizationSet(
                Casting::cast_py_to_vector<double>(
                    polarization,
                    "radian"
                )
            );
        }
        catch (const std::exception&) {
        }

        if (py::isinstance<py::list>(polarization) || py::isinstance<py::tuple>(polarization)) {
            py::sequence sequence = py::reinterpret_borrow<py::sequence>(polarization);

            if (py::len(sequence) == 0) {
                throw std::invalid_argument(
                    "Invalid type for polarization. Empty sequences are not allowed."
                );
            }

            bool all_are_polarization_state = true;

            for (size_t index = 0; index < py::len(sequence); ++index) {
                py::object item = py::reinterpret_borrow<py::object>(sequence[index]);

                if (!py::isinstance<PolarizationState>(item)) {
                    all_are_polarization_state = false;
                    break;
                }
            }

            if (all_are_polarization_state) {
                std::vector<PolarizationState> state_vector;
                state_vector.reserve(py::len(sequence));

                for (size_t index = 0; index < py::len(sequence); ++index) {
                    py::object item = py::reinterpret_borrow<py::object>(sequence[index]);
                    state_vector.push_back(item.cast<PolarizationState>());
                }

                return PolarizationSet(state_vector);
            }
        }

        throw std::invalid_argument(
            "Invalid type for polarization. Expected None, a scalar angle, a sequence of angles, "
            "a PolarizationState, or a sequence of PolarizationState."
        );
    }




    inline PolarizationState cast_py_to_polarization_state(
        const py::object& polarization
    )
    {
        if (py::isinstance<PolarizationState>(polarization)) {
            return polarization.cast<PolarizationState>();
        }

        try {
            const double angle_radian = Casting::cast_py_to_vector<double>(
                polarization,
                "radian"
            ).at(0);

            return PolarizationState(angle_radian);
        }
        catch (const std::exception&) {
            throw std::invalid_argument(
                "Invalid type for polarization. Expected a scalar angle, a Pint scalar quantity, or a PolarizationState."
            );
        }
    }










}