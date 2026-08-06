#include "labeled_array.h"
#include <pybind11/stl.h>

#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <utility>

LabeledArray::LabeledArray(
        py::object values,
        std::vector<std::string> dims,
        py::dict coords,
        py::dict attrs,
        py::object name
    ) :
        values_(py::array::ensure(values)),
        dims_(std::move(dims)),
        coords_(std::move(coords)),
        attrs_(std::move(attrs)),
        name_(std::move(name))
    {
        if (!values_) {
            throw std::invalid_argument("values must be convertible to a NumPy array");
        }
        if (static_cast<std::size_t>(values_.ndim()) != dims_.size()) {
            throw std::invalid_argument("the number of dimension names must match values.ndim");
        }

        for (std::size_t axis = 0; axis < dims_.size(); ++axis) {
            if (dims_[axis].empty()) {
                throw std::invalid_argument("dimension names cannot be empty");
            }
            if (std::find(dims_.begin(), dims_.begin() + static_cast<long>(axis), dims_[axis]) !=
                dims_.begin() + static_cast<long>(axis)) {
                throw std::invalid_argument("dimension names must be unique");
            }

            py::str dimension(dims_[axis]);
            if (!coords_.contains(dimension)) {
                py::module_ numpy = py::module_::import("numpy");
                coords_[dimension] = numpy.attr("arange")(values_.shape(axis));
            }

            py::array coordinate = py::array::ensure(coords_[dimension]);
            if (!coordinate || coordinate.ndim() != 1 || coordinate.shape(0) != values_.shape(axis)) {
                throw std::invalid_argument("each coordinate must be one-dimensional and match its axis length");
            }
            coords_[dimension] = coordinate;
        }

        for (auto item : coords_) {
            std::string dimension = py::cast<std::string>(item.first);
            if (std::find(dims_.begin(), dims_.end(), dimension) == dims_.end()) {
                throw std::invalid_argument("coordinates may only be provided for named dimensions");
            }
        }
    }

py::array LabeledArray::to_numpy() const { return values_; }

py::array LabeledArray::as_numpy() const { return values_; }

py::object LabeledArray::as_dataframe() const {
        py::module_ numpy = py::module_::import("numpy");
        py::object dataframe_type = py::module_::import("pandas").attr("DataFrame");
        py::dict columns;

        auto add_coordinate_column = [this, &numpy, &columns](
            std::size_t coordinate_axis,
            const std::string& dimension,
            py::array source_values
        ) {
            py::array coordinate = coords_[py::str(dimension)];
            py::tuple coordinate_shape(source_values.ndim());
            for (ssize_t i = 0; i < source_values.ndim(); ++i) {
                coordinate_shape[i] = py::int_(
                    i == static_cast<ssize_t>(coordinate_axis) ? source_values.shape(i) : 1
                );
            }
            py::tuple value_shape(source_values.ndim());
            for (ssize_t i = 0; i < source_values.ndim(); ++i) {
                value_shape[i] = py::int_(source_values.shape(i));
            }
            py::object reshaped = coordinate.attr("reshape")(coordinate_shape);
            py::object broadcasted = numpy.attr("broadcast_to")(reshaped, value_shape);
            columns[py::str(dimension)] = broadcasted.attr("reshape")(-1);
        };

        bool has_measure_axis = !dims_.empty() && dims_[0] == "measure";
        if (has_measure_axis) {
            py::array first_measure = py::array::ensure(numpy.attr("take")(values_, 0, py::arg("axis") = 0));
            for (std::size_t axis = 1; axis < dims_.size(); ++axis) {
                add_coordinate_column(axis - 1, dims_[axis], first_measure);
            }
            py::array measure_names = coords_[py::str("measure")];
            for (ssize_t index = 0; index < measure_names.shape(0); ++index) {
                py::object measure = measure_names.attr("__getitem__")(index);
                py::object measure_values = numpy.attr("take")(values_, index, py::arg("axis") = 0);
                columns[measure] = measure_values.attr("reshape")(-1);
            }
        } else {
            for (std::size_t axis = 0; axis < dims_.size(); ++axis) {
                add_coordinate_column(axis, dims_[axis], values_);
            }
            std::string value_name = name_.is_none() ? "value" : py::cast<std::string>(name_);
            columns[py::str(value_name)] = values_.attr("reshape")(-1);
        }

        py::object dataframe = dataframe_type(columns);
        py::dict dataframe_attrs = dataframe.attr("attrs").cast<py::dict>();
        for (auto item : attrs_) {
            dataframe_attrs[item.first] = item.second;
        }
        if (attrs_.contains(py::str("coordinate_units"))) {
            py::dict units;
            if (attrs_.contains(py::str("units"))) {
                units = attrs_.attr("__getitem__")(py::str("units")).cast<py::dict>();
            }
            py::dict coordinate_units = attrs_[py::str("coordinate_units")].cast<py::dict>();
            for (auto item : coordinate_units) {
                units[item.first] = item.second;
            }
            dataframe_attrs[py::str("units")] = units;
        }
        return dataframe;
    }

LabeledArray LabeledArray::isel(py::dict indexers) const {
        py::module_ numpy = py::module_::import("numpy");
        py::object selected_values = values_;
        std::vector<std::string> new_dims;
        py::dict new_coords;
        ssize_t current_axis = 0;

        for (std::size_t axis = 0; axis < dims_.size(); ++axis) {
            py::str dimension(dims_[axis]);
            bool provided = indexers.contains(dimension);
            if (!provided) {
                new_dims.push_back(dims_[axis]);
                new_coords[dimension] = coords_[dimension];
                ++current_axis;
                continue;
            }
            py::object selector = indexers[dimension];

            py::object selected_coordinate;
            if (py::isinstance<py::int_>(selector)) {
                ssize_t index = selector.cast<ssize_t>();
                selected_values = numpy.attr("take")(selected_values, index, py::arg("axis") = current_axis);
                selected_coordinate = numpy.attr("take")(coords_[dimension], index);
                continue;
            }
            if (!py::isinstance<py::slice>(selector)) {
                throw std::invalid_argument("isel currently accepts only integer and slice indexers");
            }

            py::tuple selectors(selected_values.attr("ndim").cast<std::size_t>());
            for (ssize_t i = 0; i < static_cast<ssize_t>(selectors.size()); ++i) {
                selectors[i] = py::slice(py::none(), py::none(), py::none());
            }
            selectors[current_axis] = selector;
            selected_values = selected_values.attr("__getitem__")(selectors);
            selected_coordinate = coords_[dimension].attr("__getitem__")(selector);
            new_dims.push_back(dims_[axis]);
            new_coords[dimension] = selected_coordinate;
            ++current_axis;
        }

        return LabeledArray(selected_values, new_dims, new_coords, attrs_, name_);
    }

py::object LabeledArray::mean(py::object dimension) const {
        py::module_ numpy = py::module_::import("numpy");
        if (dimension.is_none()) {
            return numpy.attr("mean")(values_);
        }

        std::string dim = py::cast<std::string>(dimension);
        auto iterator = std::find(dims_.begin(), dims_.end(), dim);
        if (iterator == dims_.end()) {
            throw std::invalid_argument("unknown dimension: " + dim);
        }
        std::size_t axis = static_cast<std::size_t>(std::distance(dims_.begin(), iterator));
        py::object reduced = numpy.attr("mean")(values_, py::arg("axis") = axis);
        std::vector<std::string> new_dims;
        py::dict new_coords;
        for (std::size_t i = 0; i < dims_.size(); ++i) {
            if (i != axis) {
                new_dims.push_back(dims_[i]);
                new_coords[py::str(dims_[i])] = coords_[py::str(dims_[i])];
            }
        }
        return py::cast(LabeledArray(reduced, new_dims, new_coords, attrs_, name_));
    }

py::object LabeledArray::plot(
        py::object x,
        py::object y,
        py::object ax,
        py::object cmap,
        py::object std_dimension,
        double alpha,
        bool show_legend,
        bool show
    ) const {
        if (values_.ndim() < 1) {
            throw std::invalid_argument("plot requires at least one parameter dimension");
        }

        py::module_ pyplot = py::module_::import("matplotlib.pyplot");
        if (ax.is_none()) {
            py::object figure_and_axes = pyplot.attr("subplots")();
            ax = figure_and_axes.attr("__getitem__")(1);
        }

        auto dimension_index = [this](py::object dimension) -> std::size_t {
            std::string name = py::cast<std::string>(dimension);
            auto iterator = std::find(dims_.begin(), dims_.end(), name);
            if (iterator == dims_.end()) {
                throw std::invalid_argument("unknown dimension: " + name);
            }
            return static_cast<std::size_t>(std::distance(dims_.begin(), iterator));
        };

        auto is_measure_name = [this](py::object candidate) {
            if (candidate.is_none() || !py::isinstance<py::str>(candidate)) {
                return false;
            }
            std::string name = py::cast<std::string>(candidate);
            if (!name_.is_none() && py::cast<std::string>(name_) == name) {
                return true;
            }
            if (attrs_.contains(py::str("measures"))) {
                py::iterable measures = py::reinterpret_borrow<py::iterable>(attrs_[py::str("measures")]);
                for (py::handle item : measures) {
                    if (py::cast<std::string>(item) == name) {
                        return true;
                    }
                }
            }
            return false;
        };

        auto display_label = [](py::object value) {
            std::string label = py::cast<std::string>(value);
            std::replace(label.begin(), label.end(), '_', ' ');
            return py::str(label);
        };

        auto coordinate_label = [](py::object value) {
            try {
                std::ostringstream stream;
                stream << std::setprecision(4) << value.cast<double>();
                return stream.str();
            } catch (const py::cast_error&) {
                return py::cast<std::string>(py::str(value));
            }
        };

        py::object y_measure = y.is_none() ? name_ : y;
        if (!is_measure_name(y_measure)) {
            throw std::invalid_argument("y must name a computed measure");
        }

        py::object x_dimension = x.is_none() ? py::str(dims_[0]) : x;
        std::size_t x_axis = dimension_index(x_dimension);
        py::module_ numpy = py::module_::import("numpy");
        py::list artists;
        py::array x_coordinate = coords_[py::str(x_dimension)];
        ax.attr("set_xlabel")(display_label(x_dimension));
        ax.attr("set_ylabel")(display_label(y_measure));

        bool has_measure_axis = !dims_.empty() && dims_[0] == "measure";
        py::object plot_values = values_;
        py::object spread_values = py::none();
        std::vector<std::string> plot_dims = dims_;
        std::size_t plot_x_axis = x_axis;
        std::size_t plot_std_axis = 0;

        if (has_measure_axis) {
            if (x_axis == 0) {
                throw std::invalid_argument("x must name a parameter dimension, not measure");
            }

            py::array measure_names = py::array::ensure(coords_[py::str("measure")]);
            ssize_t measure_index = -1;
            for (ssize_t index = 0; index < measure_names.shape(0); ++index) {
                if (py::cast<std::string>(measure_names.attr("__getitem__")(index)) ==
                    py::cast<std::string>(y_measure)) {
                    measure_index = index;
                    break;
                }
            }
            if (measure_index < 0) {
                throw std::invalid_argument("unknown computed measure: " + py::cast<std::string>(y_measure));
            }

            plot_values = numpy.attr("take")(values_, measure_index, py::arg("axis") = 0);
            plot_dims.erase(plot_dims.begin());
            --plot_x_axis;
        }

        if (!std_dimension.is_none()) {
            std::size_t std_axis = dimension_index(std_dimension);
            if (has_measure_axis) {
                if (std_axis == 0) {
                    throw std::invalid_argument("std must name a parameter dimension, not measure");
                }
                --std_axis;
            }
            plot_std_axis = std_axis;
            if (plot_std_axis == plot_x_axis) {
                throw std::invalid_argument("std must name a parameter different from x");
            }

            py::tuple plot_shape = plot_values.attr("shape").cast<py::tuple>();
            ssize_t std_size = plot_shape[plot_std_axis].cast<ssize_t>();
            py::object spread_source = plot_values;
            plot_values = numpy.attr("mean")(plot_values, py::arg("axis") = plot_std_axis);
            if (std_size > 1) {
                spread_values = numpy.attr("std")(
                    spread_source, py::arg("axis") = plot_std_axis, py::arg("ddof") = 1
                );
            } else {
                spread_values = numpy.attr("zeros_like")(plot_values);
            }
            plot_dims.erase(plot_dims.begin() + static_cast<long>(plot_std_axis));
            if (plot_std_axis < plot_x_axis) {
                --plot_x_axis;
            }
        }

        py::array moved_values = py::array::ensure(
            numpy.attr("moveaxis")(plot_values, plot_x_axis, 0)
        );
        py::tuple matrix_shape(2);
        matrix_shape[0] = py::int_(moved_values.shape(0));
        matrix_shape[1] = py::int_(-1);
        py::array value_matrix = py::array::ensure(moved_values.attr("reshape")(matrix_shape));

        py::array spread_matrix;
        if (!spread_values.is_none()) {
            py::array moved_spread = py::array::ensure(
                numpy.attr("moveaxis")(spread_values, plot_x_axis, 0)
            );
            spread_matrix = py::array::ensure(moved_spread.attr("reshape")(matrix_shape));
        }

        std::vector<std::size_t> other_axes;
        for (std::size_t axis = 0; axis < plot_dims.size(); ++axis) {
            if (axis != plot_x_axis) {
                other_axes.push_back(axis);
            }
        }

        for (ssize_t line_index = 0; line_index < value_matrix.shape(1); ++line_index) {
            py::object series = numpy.attr("take")(value_matrix, line_index, py::arg("axis") = 1);
            std::string label = py::cast<std::string>(display_label(y_measure));
            ssize_t remaining = line_index;
            for (ssize_t position = static_cast<ssize_t>(other_axes.size()) - 1; position >= 0; --position) {
                std::size_t axis = other_axes[static_cast<std::size_t>(position)];
                ssize_t coordinate_index = remaining % moved_values.shape(static_cast<ssize_t>(position) + 1);
                remaining /= moved_values.shape(static_cast<ssize_t>(position) + 1);
                py::object coordinate = coords_[py::str(plot_dims[axis])].attr("__getitem__")(coordinate_index);
                if (!label.empty()) {
                    label += " | ";
                }
                label += py::cast<std::string>(display_label(py::str(plot_dims[axis]))) + "=" +
                    coordinate_label(coordinate);
            }
            py::object plotted = ax.attr("plot")(
                x_coordinate, series, py::arg("label") = py::str(label)
            );
            artists.append(plotted.attr("__getitem__")(0));

            if (!spread_values.is_none()) {
                py::object spread = numpy.attr("take")(spread_matrix, line_index, py::arg("axis") = 1);
                py::object band = ax.attr("fill_between")(
                    x_coordinate,
                    series.attr("__sub__")(spread),
                    series.attr("__add__")(spread),
                    py::arg("alpha") = alpha
                );
                artists.append(band);
            }
        }
        if (show_legend) {
            ax.attr("legend")();
        }
        if (show) {
            pyplot.attr("show")();
        }
        return artists;
    }

py::tuple LabeledArray::dims() const {
        py::tuple result(dims_.size());
        for (std::size_t i = 0; i < dims_.size(); ++i) {
            result[i] = py::str(dims_[i]);
        }
        return result;
    }

py::tuple LabeledArray::shape() const {
        py::tuple result(values_.ndim());
        for (ssize_t i = 0; i < values_.ndim(); ++i) {
            result[i] = py::int_(values_.shape(i));
        }
        return result;
    }

std::string LabeledArray::repr() const {
        std::ostringstream stream;
        stream << "LabeledArray(shape=" << py::str(shape()).cast<std::string>()
               << ", dims=" << py::str(dims()).cast<std::string>() << ")";
        return stream.str();
    }
