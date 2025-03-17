# cpp/sets/CMakeLists.txt

set(NAME "experiment")

# Create a shared library for functionality.
add_library("cpp_${NAME}" SHARED "${NAME}.cpp")

# Add numpy_interface.cpp as a private source so that it’s compiled only for this target.
target_sources("cpp_${NAME}" PRIVATE "../utils/numpy_interface.h")

target_link_libraries("cpp_${NAME}" PUBLIC pybind11::module OpenMP::OpenMP_CXX cpp_cylinder cpp_sphere cpp_coreshell cpp_detector cpp_source cpp_full_mesh)


if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
    target_link_options("cpp_${NAME}" PRIVATE -lomp)
endif()

if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options("cpp_${NAME}" PRIVATE -fopenmp)
    target_link_options("cpp_${NAME}" PRIVATE -static -fopenmp -Wl,--whole-archive -lgomp -Wl,--no-whole-archive)
endif()

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "cpp_${NAME}" )



