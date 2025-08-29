# cpp/coordinates/CMakeLists.txt

set(NAME "coordinates")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

install(
    TARGETS ${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)
