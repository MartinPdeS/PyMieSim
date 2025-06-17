# cpp/full_mesh/CMakeLists.txt

set(NAME "full_mesh")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")