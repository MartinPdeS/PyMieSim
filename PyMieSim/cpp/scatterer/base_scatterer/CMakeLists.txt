# cpp/sets/CMakeLists.txt

set(NAME "base_scatterer")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")


target_link_libraries("${NAME}" PUBLIC pybind11::module OpenMP::OpenMP_CXX source full_mesh coordinates)

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xclang -fopenmp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")
    target_link_options("${NAME}" PRIVATE -lomp)
endif()

if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options("${NAME}" PRIVATE -fopenmp)
    target_link_options("${NAME}" PRIVATE -static -fopenmp -Wl,--whole-archive -lgomp -Wl,--no-whole-archive)
endif()

install(
    TARGETS ${NAME}
    LIBRARY DESTINATION ${PYMIESIM_BIN_DIR}
    ARCHIVE DESTINATION ${PYMIESIM_BIN_DIR}
    RUNTIME DESTINATION ${PYMIESIM_BIN_DIR}
)
