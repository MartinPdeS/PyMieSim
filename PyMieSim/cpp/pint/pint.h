#pragma once

#include <mutex>                // for std::once_flag, std::call_once
#include <stdexcept>            // for std::runtime_error
#include <vector>               // for std::vector
#include <string>               // for std::string
#include <pybind11/pybind11.h>  // for py::handle, py::object
#include <pybind11/stl.h>       // for py::cast with STL containers
#include <pybind11/numpy.h>     // for py::array

namespace py = pybind11;

class UnitRegistrySingleton {
public:
    static UnitRegistrySingleton& instance(); // declaration only

    UnitRegistrySingleton(const UnitRegistrySingleton&) = delete;
    UnitRegistrySingleton& operator=(const UnitRegistrySingleton&) = delete;

    void set_ureg(py::object ureg_object) {
        py::gil_scoped_acquire gil;

        std::call_once(initialization_flag, [&]() {
            ureg = std::move(ureg_object);
            if (ureg.is_none()) {
                throw std::runtime_error("UnitRegistrySingleton.set_ureg: ureg is None.");
            }
        });
    }

    py::object get_ureg() const {
        py::gil_scoped_acquire gil;

        if (ureg.is_none()) {
            throw std::runtime_error("UnitRegistrySingleton.get_ureg: ureg not initialized.");
        }
        return ureg;
    }

    bool is_initialized() const {
        py::gil_scoped_acquire gil;
        return !ureg.is_none();
    }

private:
    UnitRegistrySingleton() = default;
    mutable std::once_flag initialization_flag;
    py::object ureg = py::none();
};

/*
Function to retrieve the shared pint UnitRegistry singleton.
@returns A pint UnitRegistry object.
*/
py::object get_shared_ureg();

/*
Function to convert a pint.Quantity length to double in meters.
@param value A pint.Quantity object representing a length.
@returns double value in meters.
*/
double to_meters_strict(py::handle value);

/*
Function to retrieve the pint UnitRegistry from a pint.Quantity object.
@param value A pint.Quantity object.
@returns A pint UnitRegistry object associated with the quantity.
*/
py::object registry_from_quantity(py::handle value);

/*
Function to create a pint.Quantity object in meters from a double value.
@param ureg A pint UnitRegistry object.
@param meters_value A double value representing length in meters.
@returns A pint.Quantity object representing the length in meters.
*/
py::object meters_quantity_with_ureg(const py::object& ureg, double meters_value);

/*
Function to convert an iterable of pint.Quantity objects to a std::vector<double> in specified units.
@param values An iterable of pint.Quantity objects.
@param unit A string representing the target unit (e.g., "meter").
@returns std::vector<double> with the magnitudes in the specified unit.
*/
std::vector<double> to_vector_units(py::handle values, const std::string& unit);

/*
Function to convert a pint.Quantity scalar to double in meters.
@param quantity A pint.Quantity object representing a length.
@returns double value in meters.
*/
double quantity_scalar_to_meters(py::object quantity);

/*
Function to convert a 1D pint.Quantity array to a std::vector<double> in meters.
@param quantity A pint.Quantity object representing a 1D array of lengths.
@returns std::vector<double> with values in meters.
*/
std::vector<double> quantity_1d_to_meters_vector(py::object quantity);

/*
Function to convert an array-like object to a std::vector<double>.
@param values An array-like object (e.g., list, tuple, numpy array) containing numeric values.
@returns std::vector<double> with the numeric values.
*/
std::vector<double> array_like_1d_to_double_vector(py::object values);