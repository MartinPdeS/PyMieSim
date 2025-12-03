# cpp/fibonacci/CMakeLists.txt

set(NAME "fibonacci")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")
target_link_libraries("${NAME}" PRIVATE coordinates)

install(
    TARGETS ${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)

