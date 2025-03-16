# cpp/cylinder/CMakeLists.txt

set(NAME "cylinder")


# Create a shared library for functionality.
add_library("cpp_${NAME}" SHARED "${NAME}.cpp")

# Add numpy_interface.cpp as a private source so that it’s compiled only for this target.
target_sources("cpp_${NAME}" PRIVATE "../utils/numpy_interface.h")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module cpp_source cpp_full_mesh ZBessel _ZBessel)

target_include_directories("cpp_${NAME}" PUBLIC ${CMAKE_SOURCE_DIR}/cpp)

set_target_properties("cpp_${NAME}" PROPERTIES
    CXX_VISIBILITY_PRESET default
    VISIBILITY_INLINES_HIDDEN OFF
)

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "cpp_${NAME}" OpenMP::OpenMP_CXX)



if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    target_link_options("interface_${NAME}" PRIVATE -lomp)
endif()
