# cpp/coordinates/CMakeLists.txt
set(NAME "coordinates")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS ${NAME})
