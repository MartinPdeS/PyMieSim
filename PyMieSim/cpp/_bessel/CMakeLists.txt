# cpp/_bessel/CMakeLists.txt

# Build _ZBessel library
add_library(_ZBessel STATIC "bessel_subroutine.cpp")
target_link_options(_ZBessel PUBLIC ${STATIC_LINK_OPTIONS})