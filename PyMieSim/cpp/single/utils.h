#include <pybind11/pybind11.h>
#include <polarization/polarization.h>

namespace py = pybind11;

inline PolarizationState get_polarization_filter_state(const py::object& polarization_filter) {
    PolarizationState polarization_filter_state;

    if (polarization_filter.is_none()) {
        polarization_filter_state = PolarizationState();
    }
    else if (py::isinstance<PolarizationState>(polarization_filter)) {
        polarization_filter_state = polarization_filter.cast<PolarizationState>();
    }
    else {
        try {
            const double angle = polarization_filter.attr("to")("radian").attr("magnitude").cast<double>();

            polarization_filter_state = PolarizationState(angle);
        }
        catch (const std::exception&) {
            throw std::runtime_error(
                "Invalid type for polarization_filter. Expected None, a scalar angle, or a PolarizationState."
            );
        }
    }

    return polarization_filter_state;
}