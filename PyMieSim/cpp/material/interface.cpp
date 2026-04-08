#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include "./base.h"
#include "./constant.h"
#include "./tabulated.h"
#include "./sellmeier.h"
#include <pint/pint.h>

#include <cstdio>
#include <memory>
#include <string>
#include <vector>

namespace py = pybind11;


PYBIND11_MODULE(material, module) {
    py::object ureg = get_shared_ureg();

    module.def(
        "print_available",
        []() {
            py::module_ pyoptik = py::module_::import("PyOptik");
            py::object Material = pyoptik.attr("Material");

            Material.attr("use_tabulated") = py::bool_(true);
            Material.attr("use_sellmeier") = py::bool_(true);

            Material.attr("print_available")();
        },
        R"pbdoc(
            Prints the available materials in the PyOptik library.

            This function retrieves and prints the list of available materials from the PyOptik library. It sets the appropriate flags to indicate that both tabulated and Sellmeier materials should be included in the output.
        )pbdoc"
    );

    py::class_<BaseMaterial, std::shared_ptr<BaseMaterial>>(module, "BaseMaterial")
        .def(
            "initialize",
            [](BaseMaterial& self, const py::object& wavelength) {
                return self.initialize(
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>()
                );
            },
            py::arg("wavelength"),
            R"pbdoc(
                Initializes the material with a specific wavelength.

                This method initializes the material by computing its refractive index at a given wavelength. The wavelength is provided as an argument, and the method converts it to meters before performing the initialization. After calling this method, the material's refractive index will be available for retrieval.

                Parameters
                ----------
                wavelength : Quantity
                    The wavelength at which to initialize the material, provided as a quantity with appropriate units (e.g., nanometers, micrometers).
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const BaseMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Retrieves the refractive index for a specific wavelength.

                This method allows for querying the refractive index of the material at a specific wavelength. The wavelength is provided as an argument, and the method returns the corresponding refractive index value. If the material has not been initialized at the specified wavelength, an error will be raised.

                Parameters
                ----------
                wavelength : Quantity
                    The wavelength at which to retrieve the refractive index, provided as a quantity with appropriate units (e.g., nanometers, micrometers).
            )pbdoc"
        )
        ;

    py::class_<ConstantMaterial, BaseMaterial, std::shared_ptr<ConstantMaterial>>(module, "ConstantMaterial")
        .def(
            py::init([ureg](const py::object& refractive_index) {
                const complex128 refractive_index_value = refractive_index.cast<complex128>();

                return std::make_shared<ConstantMaterial>(refractive_index_value);
            }),
            py::arg("refractive_index"),
            R"pbdoc(
                Initializes a ConstantMaterial instance with a specified refractive index.

                This constructor allows for the creation of a ConstantMaterial instance by providing a constant refractive index value. The refractive index is expected to be provided as a complex number, which can be cast from a Python object. Once initialized, the material will have a constant refractive index that does not depend on wavelength.

                Parameters
                ----------
                refractive_index : complex
                    The constant refractive index for the material, provided as a complex number.
            )pbdoc"
        )
        .def_property_readonly(
            "refractive_index",
            [ureg](const ConstantMaterial& self) {
                py::object magnitude = py::cast(self.constant_refractive_index);
                return magnitude;
            },
            R"pbdoc(
                The constant refractive index of the material.

                This property provides access to the constant refractive index of the ConstantMaterial instance. The refractive index is returned as a complex number, which can be used for calculations and analysis involving the material's optical properties.
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const ConstantMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Retrieves the refractive index for a specific wavelength.

                This method allows for querying the refractive index of the ConstantMaterial instance at a specific wavelength. Since the material has a constant refractive index, the provided wavelength does not affect the returned value. However, this method is included for consistency with other material types and to allow for potential future extensions where the refractive index might depend on wavelength.

                Parameters
                ----------
                wavelength : Quantity
                    The wavelength at which to retrieve the refractive index, provided as a quantity with appropriate units (e.g., nanometers, micrometers). This parameter is currently ignored since the refractive index is constant, but it is included for interface consistency.
            )pbdoc"
        )
        .def(
            "__repr__",
            [](const ConstantMaterial& self) {
                char buffer[128];

                if (self.constant_refractive_index.imag() >= 0.0) {
                    std::snprintf(
                        buffer,
                        sizeof(buffer),
                        "%g + i%g",
                        self.constant_refractive_index.real(),
                        self.constant_refractive_index.imag()
                    );
                }
                else {
                    std::snprintf(
                        buffer,
                        sizeof(buffer),
                        "%g - i%g",
                        self.constant_refractive_index.real(),
                        -self.constant_refractive_index.imag()
                    );
                }

                return std::string(buffer);
            },
            R"pbdoc(
                Returns a string representation of the ConstantMaterial instance.

                This method provides a string representation of the ConstantMaterial instance, which includes the real and imaginary parts of the constant refractive index. The format is "real + i*imaginary" or "real - i*imaginary" depending on the sign of the imaginary part. This can be useful for debugging and display purposes when working with instances of ConstantMaterial.
            )pbdoc"
        )
        ;

    py::class_<TabulatedMaterial, BaseMaterial, std::shared_ptr<TabulatedMaterial>>(
            module,
            "TabulatedMaterial",
            R"pbdoc(
                TabulatedMaterial represents a material defined by tabulated refractive index data.

                This class allows for the definition of a material whose refractive index is defined by tabulated data, which consists of specific wavelength values and their corresponding refractive indices. The class provides methods to initialize the material, retrieve the refractive index for specific wavelengths, and access the material's properties such as its name, the wavelengths for which data is available, and whether extrapolation is allowed when querying outside of the defined wavelength range.

            )pbdoc"
        )
        .def(
            py::init(
                [](const py::str& material_name) {
                    py::module_ pyoptik = py::module_::import("PyOptik");
                    py::object Material = pyoptik.attr("Material");

                    Material.attr("use_tabulated") = py::bool_(true);
                    Material.attr("use_sellmeier") = py::bool_(false);

                    py::object material = Material.attr(material_name);

                    std::vector<double> n_values = material.attr("n_values").cast<std::vector<double>>();
                    std::vector<double> k_values = material.attr("k_values").cast<std::vector<double>>();

                    std::vector<complex128> refractive_index_values;
                    refractive_index_values.reserve(n_values.size());

                    for (std::size_t index = 0; index < n_values.size(); ++index) {
                        refractive_index_values.emplace_back(n_values[index], k_values[index]);
                    }

                    return std::make_shared<TabulatedMaterial>(
                        material_name.cast<std::string>(),
                        material.attr("wavelength").attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                        refractive_index_values,
                        true
                    );
                }
            ),
            py::arg("material_name"),
            R"pbdoc(
                Initializes a TabulatedMaterial instance based on a material name.

                This constructor allows for initialization of a TabulatedMaterial by providing the name of the material. The constructor retrieves the necessary data from the PyOptik library to set up the TabulatedMaterial instance accordingly.

                Parameters
                ----------
                material_name : str
                    The name of the material to initialize, which should correspond to a material defined in the PyOptik library.
            )pbdoc"
        )
        .def(
            py::init(
                [ureg](
                    const std::string& name,
                    const py::object& wavelengths,
                    const py::object& refractive_indices,
                    const bool allow_extrapolation
                ) {
                    const std::vector<double> wavelength_values =
                        wavelengths.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                    const std::vector<complex128> refractive_index_values = refractive_indices.cast<std::vector<complex128>>();

                    return std::make_shared<TabulatedMaterial>(
                        name,
                        wavelength_values,
                        refractive_index_values,
                        allow_extrapolation
                    );
                }
            ),
            py::arg("name"),
            py::arg("wavelengths"),
            py::arg("refractive_indices"),
            py::arg("allow_extrapolation") = false,
            R"pbdoc(
                Initializes a TabulatedMaterial instance with specified parameters.

                This constructor allows for direct initialization of a TabulatedMaterial by providing the necessary parameters, including the name, wavelengths, refractive indices, and whether extrapolation is allowed.

                Parameters
                ----------
                name : str
                    The name of the material.
                wavelengths : List[float]
                    The wavelengths (in meters) for which the refractive index data is provided.
                refractive_indices : List[complex]
                    The refractive indices corresponding to the provided wavelengths.
                allow_extrapolation : bool, optional
                    Indicates whether extrapolation is allowed when computing the refractive index outside of the defined wavelength range. Default is False.
            )pbdoc"
        )
        .def_readonly(
            "name",
            &TabulatedMaterial::name,
            R"pbdoc(
                The name of the material.

                This attribute holds the name of the TabulatedMaterial instance, which can be used for identification and representation purposes.
            )pbdoc"
        )
        .def_property_readonly(
            "wavelengths",
            [ureg](const TabulatedMaterial& self) {
                py::object magnitude = py::cast(self.wavelengths);
                return magnitude * ureg.attr("meter");
            },
            R"pbdoc(
                The wavelengths for which the refractive index data is provided.

                This property provides access to the wavelengths (in meters) for which the refractive index data of the TabulatedMaterial instance is defined. The wavelengths are returned as a list of values with appropriate units.
            )pbdoc"
        )
        .def_property_readonly(
            "refractive_indices",
            [ureg](const TabulatedMaterial& self) {
                py::object magnitude = py::cast(self.refractive_indices);
                return magnitude;
            },
            R"pbdoc(
                The refractive indices corresponding to the provided wavelengths.

                This property provides access to the refractive index values that correspond to the wavelengths defined in the `wavelengths` property. The refractive indices are returned as a list of complex values, which can be used for analysis and interpolation of the material's optical properties.
            )pbdoc"
        )
        .def_readonly(
            "allow_extrapolation",
            &TabulatedMaterial::allow_extrapolation,
            R"pbdoc(
                Indicates whether extrapolation is allowed when computing the refractive index outside of the defined wavelength range.

                This attribute indicates whether the TabulatedMaterial instance allows for extrapolation when querying the refractive index for wavelengths that fall outside of the defined range of wavelengths. If `allow_extrapolation` is set to `True`, the material will provide an extrapolated refractive index value based on the available data; if set to `False`, an error will be raised when attempting to query outside of the defined wavelength range.
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const TabulatedMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Retrieves the refractive index for a specific wavelength.

                This method allows for querying the refractive index of the TabulatedMaterial instance at a specific wavelength. The wavelength is provided as an argument, and the method returns the corresponding refractive index value. If the wavelength falls outside of the defined range and `allow_extrapolation` is set to `False`, an error will be raised.
            )pbdoc"
        )
        .def(
            "__repr__",
            [](const TabulatedMaterial& self) {
                return self.name;
            },
            R"pbdoc(
                Returns a string representation of the TabulatedMaterial instance.

                This method provides a string representation of the TabulatedMaterial instance, which in this case is simply the name of the material. This can be useful for debugging and display purposes when working with instances of TabulatedMaterial.
            )pbdoc"
        )
        ;

    py::class_<SellmeierMaterial, BaseMaterial, std::shared_ptr<SellmeierMaterial>>(
            module,
            "SellmeierMaterial",
            R"pbdoc(
                SellmeierMaterial represents a material defined by the Sellmeier equation.

                This class allows for the definition of a material whose refractive index is computed based on the Sellmeier equation, which is a common empirical relationship used to describe the wavelength dependence of the refractive index of transparent materials. The class provides methods to initialize the material and to get the refractive index for specific wavelengths.
            )pbdoc"
        )
        .def(
            py::init([](const py::str& material_name) {
                py::module_ pyoptik = py::module_::import("PyOptik");
                py::object Material = pyoptik.attr("Material");

                Material.attr("use_tabulated") = py::bool_(false);
                Material.attr("use_sellmeier") = py::bool_(true);

                py::object material = Material.attr(material_name);

                return std::make_shared<SellmeierMaterial>(
                    material_name.cast<std::string>(),
                    material.attr("coefficients").cast<std::vector<double>>(),
                    material.attr("formula_type").cast<std::size_t>(),
                    material.attr("wavelength_bound").attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                    true
                );
            }),
            py::arg("material_name"),
            R"pbdoc(
                Initializes a SellmeierMaterial instance based on a material name.

                This constructor allows for initialization of a SellmeierMaterial by providing the name of the material. The constructor retrieves the necessary data from the PyOptik library to set up the SellmeierMaterial instance accordingly.

                Parameters
                ----------
                material_name : str
                    The name of the material to initialize, which should correspond to a material defined in the PyOptik library.
            )pbdoc"
        )
        .def(
            py::init(
                [ureg](
                    const std::string& name,
                    const std::vector<double>& coefficients,
                    const std::size_t formula_type,
                    const py::object& wavelength_bound,
                    const bool allow_extrapolation
                ) {
                    const std::vector<double> wavelength_bound_values =
                        wavelength_bound.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                    return std::make_shared<SellmeierMaterial>(
                        name,
                        coefficients,
                        formula_type,
                        wavelength_bound_values,
                        allow_extrapolation
                    );
                }
            ),
            py::arg("name"),
            py::arg("coefficients"),
            py::arg("formula_type"),
            py::arg("wavelength_bound") = py::cast(std::vector<double>{}) * ureg.attr("meter"),
            py::arg("allow_extrapolation") = false,
            R"pbdoc(
                Initializes a SellmeierMaterial instance with specified parameters.

                This constructor allows for direct initialization of a SellmeierMaterial by providing the necessary parameters, including the name, coefficients, formula type, wavelength bounds, and whether extrapolation is allowed.

                Parameters
                ----------
                name : str
                    The name of the material.
                coefficients : List[float]
                    The coefficients used in the Sellmeier equation to compute the refractive index.
                formula_type : int
                    The type of Sellmeier formula to use for computing the refractive index.
                wavelength_bound : List[float], optional
                    The wavelength bounds (in meters) for which the Sellmeier equation is valid. If not provided, no bounds are applied.
                allow_extrapolation : bool, optional
                    Indicates whether extrapolation is allowed when computing the refractive index outside of the defined wavelength bounds. Default is False.
            )pbdoc"
        )
        .def_readonly(
            "name",
            &SellmeierMaterial::name,
            R"pbdoc(
                The name of the material.

                This attribute provides access to the name of the SellmeierMaterial instance, which can be used for identification and representation purposes.
            )pbdoc"
        )
        .def_readonly(
            "coefficients",
            &SellmeierMaterial::coefficients,
            R"pbdoc(
                The coefficients used in the Sellmeier equation to compute the refractive index.

                This attribute provides access to the coefficients of the SellmeierMaterial instance, which are used in the Sellmeier equation to calculate the refractive index based on the wavelength.
            )pbdoc"
        )
        .def_readonly(
            "formula_type",
            &SellmeierMaterial::formula_type,
            R"pbdoc(
                The type of Sellmeier formula to use for computing the refractive index.

                This attribute provides access to the formula type of the SellmeierMaterial instance, which determines the specific Sellmeier equation used to calculate the refractive index.
            )pbdoc"
        )
        .def_readonly(
            "allow_extrapolation",
            &SellmeierMaterial::allow_extrapolation,
            R"pbdoc(
                Indicates whether extrapolation is allowed when computing the refractive index outside of the defined wavelength bounds.

                This attribute provides access to the allow_extrapolation flag of the SellmeierMaterial instance, which determines whether the refractive index can be extrapolated beyond the specified wavelength bounds.
            )pbdoc"
        )
        .def_property_readonly(
            "wavelength_bound",
            [ureg](const SellmeierMaterial& self) {
                py::object magnitude = py::cast(self.wavelength_bound);
                return magnitude * ureg.attr("meter");
            },
            R"pbdoc(
                The wavelength bounds for which the Sellmeier equation is valid.

                This property provides access to the wavelength bounds of the SellmeierMaterial instance, which indicate the range of wavelengths (in meters) for which the Sellmeier equation is applicable. The bounds are returned as a quantity with units of meters.
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const SellmeierMaterial& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Retrieves the refractive index for a given wavelength.

                This method takes a wavelength as input and returns the corresponding refractive index of the material based on the Sellmeier equation. The wavelength is expected to be a quantity that can be converted to meters, and the method computes the refractive index accordingly.
            )pbdoc"
        )
        .def(
            "__repr__",
            [](const SellmeierMaterial& self) {
                return self.name;
            },
            R"pbdoc(
                Returns a string representation of the SellmeierMaterial instance.

                This method provides a string representation of the SellmeierMaterial instance, which includes the name of the material. This can be useful for debugging and visualization purposes.
            )pbdoc"
        )
        ;

    py::class_<BaseMedium, std::shared_ptr<BaseMedium>>(
            module,
            "BaseMedium",
            R"pbdoc(
                BaseMedium is an abstract base class representing a medium with a refractive index.

                This class defines the interface for materials that have a refractive index, which can be initialized based on a wavelength and retrieved for specific wavelengths. The class provides methods to initialize the medium and to get the refractive index, which must be implemented by derived classes.
            )pbdoc"
        )
        .def(
            "initialize",
            [](BaseMedium& self, const py::object& wavelength) {
                return self.initialize(
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>()
                );
            },
            py::arg("wavelength"),
            R"pbdoc(
                Initializes the medium with a specific wavelength.

                This method computes the refractive index of the medium based on the provided wavelength and sets the initialized state of the medium. The wavelength is expected to be a quantity that can be converted to meters, and the refractive index is computed accordingly.
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const BaseMedium& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Retrieves the refractive index for a given wavelength.

                This method takes a wavelength as input and returns the corresponding refractive index of the medium. The wavelength is expected to be a quantity that can be converted to meters, and the method computes the refractive index based on the initialized state of the medium.
            )pbdoc"
        )
        ;

    py::class_<ConstantMedium, BaseMedium, std::shared_ptr<ConstantMedium>>(
            module,
            "ConstantMedium",
            R"pbdoc(
                ConstantMedium represents a medium with a constant refractive index.

                This class allows for the definition of a medium with a constant refractive index, which does not depend on the wavelength. The class provides methods to initialize the medium and to retrieve the refractive index, which will always return the same value regardless of the wavelength.

            )pbdoc"
        )
        .def(
            py::init([ureg](const py::object& refractive_index) {
                const double refractive_index_value = refractive_index.cast<double>();

                return std::make_shared<ConstantMedium>(refractive_index_value);
            }),
            py::arg("refractive_index"),
            R"pbdoc(
                Initializes a ConstantMedium instance with a specified refractive index.

                This constructor allows for direct initialization of a ConstantMedium by providing the constant refractive index value.

                Parameters
                ----------
                refractive_index : float
                    The constant refractive index of the medium.
            )pbdoc"
        )
        .def_property_readonly(
            "refractive_index",
            [ureg](const ConstantMedium& self) {
                py::object magnitude = py::cast(self.constant_refractive_index);
                return magnitude;
            },
            R"pbdoc(
                The constant refractive index of the medium.

                This property provides access to the constant refractive index value of the ConstantMedium instance, which can be used for analysis and representation purposes.
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const ConstantMedium& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Retrieves the refractive index for a given wavelength.

                This method takes a wavelength as input and returns the constant refractive index of the medium. Since the refractive index is constant, the input wavelength does not affect the output value.
            )pbdoc"
        )
        .def(
            "__repr__",
            [](const ConstantMedium& self) {
                char buffer[64];
                std::snprintf(buffer, sizeof(buffer), "%g", self.constant_refractive_index);
                return std::string(buffer);
            },
            R"pbdoc(
                Returns a string representation of the ConstantMedium instance.

                This method provides a string representation of the ConstantMedium instance, which includes the constant refractive index value. This can be useful for debugging and visualization purposes.
            )pbdoc"
        )
        ;

    py::class_<TabulatedMedium, BaseMedium, std::shared_ptr<TabulatedMedium>>(
            module,
            "TabulatedMedium",
            R"pbdoc(
                TabulatedMedium represents a medium whose refractive index is defined by tabulated data.

                This class allows for the definition of a medium based on tabulated values of refractive index at specific wavelengths. The refractive index can be interpolated for wavelengths within the range of the provided data, and optionally extrapolated outside that range if allowed. The class provides methods to initialize the medium and to retrieve the refractive index for given wavelengths.
            )pbdoc"
        )
        .def(
            py::init(
                [ureg](
                    const std::string& name,
                    const py::object& wavelengths,
                    const py::object& refractive_indices,
                    const bool allow_extrapolation
                ) {
                    const std::vector<double> wavelength_values = wavelengths.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                    const std::vector<double> refractive_index_values = refractive_indices.cast<std::vector<double>>();

                    return std::make_shared<TabulatedMedium>(
                        name,
                        wavelength_values,
                        refractive_index_values,
                        allow_extrapolation
                    );
                }
            ),
            py::arg("name"),
            py::arg("wavelengths"),
            py::arg("refractive_indices"),
            py::arg("allow_extrapolation") = false,
            R"pbdoc(
                Initializes a TabulatedMedium instance with specified parameters.

                This constructor allows for direct initialization of a TabulatedMedium by providing the necessary parameters such as the name, wavelengths, refractive indices, and whether to allow extrapolation.

                Parameters
                ----------
                name : str
                    The name of the tabulated medium.
                wavelengths : List[float]
                    The list of wavelengths corresponding to the tabulated refractive index values, in meters.
                refractive_indices : List[complex]
                    The list of refractive index values corresponding to the wavelengths.
                allow_extrapolation : bool
                    Whether to allow extrapolation of the refractive index outside the range of the provided wavelengths. Default is False.
            )pbdoc"
        )
        .def_readonly(
            "name",
            &TabulatedMedium::name,
            R"pbdoc(
                The name of the tabulated medium.

                This property provides access to the name of the TabulatedMedium instance, which can be used for identification and representation purposes.
            )pbdoc"
        )
        .def_property_readonly(
            "wavelengths",
            [ureg](const TabulatedMedium& self) {
                py::object magnitude = py::cast(self.wavelengths);
                return magnitude * ureg.attr("meter");
            },
            R"pbdoc(
                The wavelengths corresponding to the tabulated refractive index values.

                This property returns a list of wavelengths for which the refractive index values are provided. The wavelengths are returned in meters, and can be used to understand the range of data available for the TabulatedMedium.
            )pbdoc"
        )
        .def_property_readonly(
            "refractive_indices",
            [ureg](const TabulatedMedium& self) {
                py::object magnitude = py::cast(self.refractive_indices);
                return magnitude;
            },
            R"pbdoc(
                The refractive index values corresponding to the wavelengths.

                This property returns a list of complex refractive index values that correspond to the wavelengths provided in the `wavelengths` property. These values can be used for interpolation and analysis of the medium's optical properties.
            )pbdoc"
        )
        .def_readonly(
            "allow_extrapolation",
            &TabulatedMedium::allow_extrapolation,
            R"pbdoc(
                Indicates whether extrapolation of the refractive index is allowed.

                This property indicates whether the TabulatedMedium instance allows for extrapolation of the refractive index values outside the range of the provided wavelengths. If set to True, the `get_refractive_index` method will return extrapolated values for wavelengths outside the data range; if False, it will raise an error for such cases.
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const TabulatedMedium& self, const py::object& wavelength) {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Retrieves the refractive index for a given wavelength.

                This method takes a wavelength as input and returns the corresponding refractive index. If the wavelength is within the range of the provided data, it will return an interpolated value. If the wavelength is outside the range and extrapolation is allowed, it will return an extrapolated value; otherwise, it will raise an error.

                Parameters
                ----------
                wavelength : float
                    The wavelength for which to retrieve the refractive index, in meters.

                Returns
                -------
                complex
                    The refractive index corresponding to the given wavelength.
            )pbdoc"
        )
        .def(
            "__repr__",
            [](const TabulatedMedium& self) {
                return self.name;
            },
            R"pbdoc(
                Returns a string representation of the TabulatedMedium.

                This method provides a string representation of the TabulatedMedium instance, which in this case is simply the name of the medium. This can be useful for debugging and logging purposes, allowing for easy identification of the material when printed or displayed.
            )pbdoc"
        )
        ;

    py::class_<SellmeierMedium, BaseMedium, std::shared_ptr<SellmeierMedium>>(
            module,
            "SellmeierMedium",
                R"pbdoc(
                    SellmeierMedium represents a medium whose refractive index is computed using a Sellmeier formula.

                    The Sellmeier formula is an empirical relationship that describes how the refractive index of a material varies with wavelength. It is commonly used in optics to model the dispersion of light in transparent materials.

                    The class provides methods to initialize the medium with specific parameters and to compute the refractive index for given wavelengths. It also includes properties to access the material's name, coefficients, formula type, and wavelength bounds.
                )pbdoc"
        )
        .def(
            py::init(
                [](const py::str& material_name) {
                    py::module_ pyoptik = py::module_::import("PyOptik");
                    py::object Material = pyoptik.attr("Material");

                    Material.attr("use_tabulated") = py::bool_(false);
                    Material.attr("use_sellmeier") = py::bool_(true);

                    py::object material = Material.attr(material_name);

                    return std::make_shared<SellmeierMedium>(
                        material_name.cast<std::string>(),
                        material.attr("coefficients").cast<std::vector<double>>(),
                        material.attr("formula_type").cast<std::size_t>(),
                        material.attr("wavelength_bound").attr("to")("meter").attr("magnitude").cast<std::vector<double>>(),
                        true
                    );
                }
            ),
            py::arg("material_name"),
            R"pbdoc(
                Initializes a SellmeierMedium instance based on a material name.

                This constructor retrieves the material properties from the PyOptik library using the provided material name. It initializes the SellmeierMedium with the corresponding coefficients, formula type, and wavelength bounds defined in PyOptik.

                Parameters
                ----------
                material_name : str
                    The name of the material to initialize the SellmeierMedium with. This name should correspond to a material defined in the PyOptik library.
            )pbdoc"
        )
        .def(
            py::init([ureg](
                const std::string& name,
                const std::vector<double>& coefficients,
                const std::size_t formula_type,
                const py::object& wavelength_bound,
                const bool allow_extrapolation
            ) {
                const std::vector<double> wavelength_bound_values =
                    wavelength_bound.attr("to")("meter").attr("magnitude").cast<std::vector<double>>();

                return std::make_shared<SellmeierMedium>(
                    name,
                    coefficients,
                    formula_type,
                    wavelength_bound_values,
                    allow_extrapolation
                );
            }),
            py::arg("name"),
            py::arg("coefficients"),
            py::arg("formula_type"),
            py::arg("wavelength_bound") = py::cast(std::vector<double>{}) * ureg.attr("meter"),
            py::arg("allow_extrapolation") = false,
            R"pbdoc(
                Initializes a SellmeierMedium instance with specified parameters.

                This constructor allows for direct initialization of a SellmeierMedium by providing the necessary parameters such as the name, coefficients, formula type, wavelength bounds, and whether to allow extrapolation.

                Parameters
                ----------
                name : str
                    The name of the Sellmeier medium.
                coefficients : List[float]
                    The coefficients used in the Sellmeier formula.
                formula_type : int
                    The type of Sellmeier formula to use (e.g., 1, 2, 5, or 6).
                wavelength_bound : Optional[List[float]]
                    The lower and upper bounds of the wavelength range for which the Sellmeier formula is valid, in meters. If not provided, there are no bounds.
                allow_extrapolation : bool
                    Whether to allow extrapolation outside the wavelength bounds. Default is False.
            )pbdoc"
        )
        .def_readonly(
            "name",
            &SellmeierMedium::name,
            R"pbdoc(
                The name of the Sellmeier medium.

                This property provides access to the name of the Sellmeier medium, which is a string identifier for the material. The name can be used for reference and identification purposes when working with multiple materials in simulations or calculations.
            )pbdoc"
        )
        .def_readonly(
            "coefficients",
            &SellmeierMedium::coefficients,
            R"pbdoc(
                The coefficients used in the Sellmeier formula.

                This property provides access to the coefficients of the Sellmeier medium, which are used in the Sellmeier formula to calculate the refractive index. The coefficients are typically provided as a list of floating-point numbers.
            )pbdoc"

        )
        .def_readonly(
            "formula_type",
            &SellmeierMedium::formula_type,
            R"pbdoc(
                The type of Sellmeier formula to use.

                This property indicates the specific type of Sellmeier formula that should be applied when calculating the refractive index. Common types include 1, 2, 5, and 6, each corresponding to a different mathematical form of the Sellmeier equation.
            )pbdoc"
        )
        .def_readonly(
            "allow_extrapolation",
            &SellmeierMedium::allow_extrapolation,
            R"pbdoc(
                Whether to allow extrapolation outside the wavelength bounds.

                This property indicates whether the SellmeierMedium should allow extrapolation of the refractive index calculation beyond the specified wavelength bounds. If set to True, the medium will provide refractive index values even for wavelengths outside the defined range, which may lead to less accurate results. If set to False, an error will be raised when attempting to calculate the refractive index for out-of-bounds wavelengths.
            )pbdoc"
        )
        .def_property_readonly(
            "wavelength_bound",
            [ureg](const SellmeierMedium& self) -> py::object {
                py::object magnitude = py::cast(self.wavelength_bound);
                return magnitude * ureg.attr("meter");
            },
            R"pbdoc(
                The wavelength bounds for the Sellmeier formula.

                This property provides access to the lower and upper bounds of the wavelength range for which the Sellmeier formula is valid. The bounds are returned as a list of two floating-point numbers, representing the minimum and maximum wavelengths in meters. If no bounds were specified during initialization, this property will return an empty list.
            )pbdoc"
        )
        .def(
            "get_refractive_index",
            [ureg](const SellmeierMedium& self, const py::object& wavelength) -> py::object {
                const double wavelength_meter =
                    wavelength.attr("to")("meter").attr("magnitude").cast<double>();

                py::object magnitude = py::cast(self.get_refractive_index(wavelength_meter));
                return magnitude;
            },
            py::arg("wavelength"),
            R"pbdoc(
                Computes the refractive index for a given wavelength.

                This method calculates the refractive index of the SellmeierMedium at a specified wavelength using the defined Sellmeier formula and coefficients. The input wavelength should be provided as a quantity with units (e.g., using Pint), and the method will handle the conversion to meters internally. The output is the computed refractive index, which may be a complex number depending on the material properties.

                Parameters
                ----------
                wavelength : Quantity
                    The wavelength at which to compute the refractive index. This should be a quantity with appropriate units (e.g., nanometers, micrometers, etc.) that can be converted to meters.
            )pbdoc"
        )
        .def(
            "__repr__",
            [](const SellmeierMedium& self) -> py::str {
                return self.name;
            },
            R"pbdoc(
                Returns a string representation of the SellmeierMedium.

                This method provides a string representation of the SellmeierMedium instance, which in this case is simply the name of the medium. This can be useful for debugging and logging purposes, allowing for easy identification of the material when printed or displayed.
            )pbdoc"
        )
        ;
}