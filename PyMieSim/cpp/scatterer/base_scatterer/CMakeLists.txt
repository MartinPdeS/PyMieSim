# cpp/sets/CMakeLists.txt

set(NAME "base_scatterer")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module source full_mesh coordinates)

install(
    TARGETS ${NAME}
    LIBRARY DESTINATION ${PYMIESIM_BIN_DIR}
    ARCHIVE DESTINATION ${PYMIESIM_BIN_DIR}
    RUNTIME DESTINATION ${PYMIESIM_BIN_DIR}
)