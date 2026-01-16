# cpp/experiment/CMakeLists.txt

set(NAME "experiment")

add_subdirectory(setup)

add_subdirectory(sets)

# Create a Python module.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PUBLIC sets setup)


set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "interface_${NAME}")