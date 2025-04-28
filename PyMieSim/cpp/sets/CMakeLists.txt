# cpp/sets/CMakeLists.txt

set(NAME "sets")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module cpp_sphere cpp_cylinder cpp_coreshell cpp_source cpp_full_mesh ZBessel)

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "cpp_${NAME}")
