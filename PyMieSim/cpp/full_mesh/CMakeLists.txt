# cpp/full_mesh/CMakeLists.txt

set(NAME "full_mesh")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "${NAME}")