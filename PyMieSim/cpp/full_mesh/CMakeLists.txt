# cpp/full_mesh/CMakeLists.txt

set(NAME "full_mesh")

# Create a shared library for functionality.
add_library("cpp_${NAME}" SHARED "${NAME}.cpp")

# Add numpy_interface.cpp as a private source so that it’s compiled only for this target.
target_sources("cpp_${NAME}" PRIVATE "../utils/numpy_interface.h")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module)

target_include_directories("cpp_${NAME}" PUBLIC ${CMAKE_SOURCE_DIR}/cpp)

set_target_properties("cpp_${NAME}" PROPERTIES
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
)
