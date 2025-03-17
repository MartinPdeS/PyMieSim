# cpp/full_mesh/CMakeLists.txt

set(NAME "full_mesh")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")

# Add numpy_interface.cpp as a private source so that it’s compiled only for this target.
target_sources("cpp_${NAME}" PRIVATE "../utils/numpy_interface.h")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module)