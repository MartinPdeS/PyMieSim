# cpp/sphere/CMakeLists.txt

set(NAME "sphere")

# Create a shared library for functionality.
add_library("cpp_${NAME}" SHARED "${NAME}.cpp")

# Add numpy_interface.cpp as a private source so that it’s compiled only for this target.
target_sources("cpp_${NAME}" PRIVATE "../utils/numpy_interface.h")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module cpp_source ZBessel PUBLIC cpp_full_mesh)

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "cpp_${NAME}")
