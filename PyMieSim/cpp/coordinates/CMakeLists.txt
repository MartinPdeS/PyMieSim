# cpp/coordinates/CMakeLists.txt

set(NAME "coordinates")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")

install(
    TARGETS cpp_${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)