# cpp/detector/CMakeLists.txt

set(NAME "detector")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp" "../utils/numpy_interface.h")

target_link_libraries("cpp_${NAME}" PUBLIC pybind11::module ZBessel PRIVATE cpp_mode_field cpp_fibonacci cpp_cylinder cpp_sphere cpp_coreshell cpp_full_mesh)

target_include_directories("cpp_${NAME}" PUBLIC ${CMAKE_SOURCE_DIR}/cpp)

set_target_properties("cpp_${NAME}" PROPERTIES
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
)

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "cpp_${NAME}")