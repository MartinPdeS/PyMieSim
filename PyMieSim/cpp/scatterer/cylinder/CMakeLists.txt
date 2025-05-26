# cpp/cylinder/CMakeLists.txt

set(NAME "cylinder")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("cpp_${NAME}" PRIVATE pybind11::module ZBessel PUBLIC cpp_source cpp_full_mesh cpp_base_scatterer)

install(
    TARGETS cpp_${NAME}
    LIBRARY DESTINATION ${PYMIESIM_BIN_DIR}
    ARCHIVE DESTINATION ${PYMIESIM_BIN_DIR}
    RUNTIME DESTINATION ${PYMIESIM_BIN_DIR}
)