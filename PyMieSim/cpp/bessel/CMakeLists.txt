# cpp/bessel/CMakeLists.txt

# Build ZBessel library
add_library(ZBessel STATIC
    "amos_iso_c_fortran_wrapper.f90"
    "machine.for"
    "zbesh.for"
)
target_link_options(ZBessel PUBLIC ${STATIC_LINK_OPTIONS})

install(
    TARGETS ZBessel
    LIBRARY DESTINATION ${PYMIESIM_BIN_DIR}
    ARCHIVE DESTINATION ${PYMIESIM_BIN_DIR}
    RUNTIME DESTINATION ${PYMIESIM_BIN_DIR}
)