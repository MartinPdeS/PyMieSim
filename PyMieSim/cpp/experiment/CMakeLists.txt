# cpp/sets/CMakeLists.txt

set(NAME "experiment")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module OpenMP::OpenMP_CXX cylinder sphere coreshell detector source full_mesh)

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
    target_link_options("${NAME}" PRIVATE -lomp)
endif()

if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options("${NAME}" PRIVATE -fopenmp)
    target_link_options("${NAME}" PRIVATE -static -fopenmp -Wl,--whole-archive -lgomp -Wl,--no-whole-archive)
endif()

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "${NAME}" )


install(
    TARGETS ${NAME} interface_${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)
