# cpp/single/CMakeLists.txt

set(NAME "single")


add_subdirectory(scatterer/base_scatterer)     # base_scatterer
add_subdirectory(scatterer/sphere)             # sphere         - interface_sphere
add_subdirectory(scatterer/cylinder)           # cylinder       - interface_cylinder
add_subdirectory(scatterer/coreshell)          # coreshell      - interface_coreshell

add_subdirectory(source)
add_subdirectory(detector)

# Create a Python module.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC detector sphere cylinder coreshell base_scatterer full_mesh mode_field fibonacci)


set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "interface_${NAME}")