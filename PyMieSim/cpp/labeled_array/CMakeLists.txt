set(NAME "labeled_array")

pybind11_add_module("interface_${NAME}" MODULE interface.cpp labeled_array.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC pybind11::module)

install(
    TARGETS "interface_${NAME}"
    LIBRARY DESTINATION "PyMieSim"
    RUNTIME DESTINATION "PyMieSim"
)
