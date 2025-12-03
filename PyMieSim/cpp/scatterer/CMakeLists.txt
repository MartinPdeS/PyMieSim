# cpp/sets/CMakeLists.txt

set(NAME "scatterer")


add_subdirectory(base_scatterer)     # base_scatterer
add_subdirectory(sphere)             # sphere         - interface_sphere
add_subdirectory(cylinder)           # cylinder       - interface_cylinder
add_subdirectory(coreshell)          # coreshell      - interface_coreshell


# Create a Python module.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC sphere cylinder coreshell base_scatterer)


set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "interface_${NAME}")
