# cpp/detector/CMakeLists.txt

set(NAME "detector")

add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module mode_field fibonacci cylinder sphere coreshell full_mesh optical_interface)
