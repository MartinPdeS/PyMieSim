#include "utils.h"
#include <pint/pint.h>

std::shared_ptr<BaseMaterial> parse_material_object(
    const py::object& material_object,
    const py::object& ureg
) {
    if (py::isinstance<BaseMaterial>(material_object)) {
        return material_object.cast<std::shared_ptr<BaseMaterial>>();
    }

    validated_no_unit(material_object, "material_object", true);

    if (py::isinstance<py::float_>(material_object) || py::isinstance<py::int_>(material_object)) {
        const complex128 refractive_index_value =
            py::cast<complex128>(material_object);

        return std::make_shared<ConstantMaterial>(refractive_index_value);
    }

    try {
        const complex128 refractive_index_value =
            py::cast<complex128>(material_object);

        return std::make_shared<ConstantMaterial>(refractive_index_value);
    }
    catch (const py::cast_error&) {
    }


    throw std::runtime_error(
        "material must be either a BaseMaterial instance or a constant complex refractive index."
    );
}

std::shared_ptr<BaseMedium> parse_medium_object(
    const py::object& medium_object,
    const py::object& ureg
) {
    if (py::isinstance<BaseMedium>(medium_object)) {
        return medium_object.cast<std::shared_ptr<BaseMedium>>();
    }

    validated_no_unit(medium_object, "medium_object", true);

    if (py::isinstance<py::float_>(medium_object) || py::isinstance<py::int_>(medium_object)) {
        const double refractive_index_value =
            py::cast<double>(medium_object);

        return std::make_shared<ConstantMedium>(refractive_index_value);
    }


    throw std::runtime_error(
        "medium must be either a BaseMedium instance or a constant real refractive index."
    );
}
