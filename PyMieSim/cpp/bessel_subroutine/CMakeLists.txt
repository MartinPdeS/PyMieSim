# cpp/bessel/CMakeLists.txt



# Build bessel_subroutine library
add_library(bessel_subroutine STATIC
    "amos_iso_c_fortran_wrapper.f90"
    "machine.for"
    "zbesh.for"
)
target_link_options(bessel_subroutine PUBLIC ${STATIC_LINK_OPTIONS})


target_compile_options(bessel_subroutine PRIVATE -Wno-maybe-uninitialized)

install(
    TARGETS bessel_subroutine
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)
