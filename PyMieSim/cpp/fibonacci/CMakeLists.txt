# cpp/fibonacci/CMakeLists.txt

set(NAME "fibonacci")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")
target_link_libraries("${NAME}" PRIVATE coordinates)

set_property(GLOBAL APPEND PROPERTY PYMIESIM_TARGETS "${NAME}")

