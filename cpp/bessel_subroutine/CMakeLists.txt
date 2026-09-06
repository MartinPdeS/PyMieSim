# cpp/bessel/CMakeLists.txt
set(NAME "bessel_subroutine")

# Build bessel_subroutine library
add_library(${NAME} STATIC
    "errors.cpp"
    "amos_iso_c_fortran_wrapper.f90"
    "machine.for"
    "zbesh.for"
)
target_link_options(${NAME} PUBLIC ${STATIC_LINK_OPTIONS})
