# cpp/coreshell/CMakeLists.txt

set(NAME "coreshell")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module bessel_subroutine PUBLIC source full_mesh base_scatterer)

install(
    TARGETS ${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)
