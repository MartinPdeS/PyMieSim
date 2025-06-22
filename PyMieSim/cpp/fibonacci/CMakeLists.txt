# cpp/fibonacci/CMakeLists.txt

set(NAME "fibonacci")

# Create a shared library for functionality.
add_library("cpp_${NAME}" STATIC "${NAME}.cpp")
target_link_libraries("cpp_${NAME}" PRIVATE cpp_coordinates)

install(
    TARGETS cpp_${NAME}
    LIBRARY DESTINATION ${LOCAL_BIN_DIR}
    ARCHIVE DESTINATION ${LOCAL_BIN_DIR}
    RUNTIME DESTINATION ${LOCAL_BIN_DIR}
)