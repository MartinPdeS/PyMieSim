# cpp/cylinder/CMakeLists.txt

set(NAME "cylinder")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PRIVATE pybind11::module bessel_subroutine PUBLIC source full_mesh base_scatterer)

install(
    TARGETS ${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)