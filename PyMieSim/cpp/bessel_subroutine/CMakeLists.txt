# cpp/bessel_subroutine/CMakeLists.txt
set(NAME "bessel_subroutine")



# Build bessel_subroutine library
add_library(${NAME} STATIC
    # "amos_iso_c_fortran_wrapper.f90"
    # "machine.for"
    # "zbesh.for"
)
target_link_options("${NAME}" PUBLIC "${STATIC_LINK_OPTIONS}")

target_compile_options(${NAME} PRIVATE -Wno-maybe-uninitialized)

install(
    TARGETS ${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)
