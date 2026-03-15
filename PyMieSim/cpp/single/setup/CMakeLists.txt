# PyMieSim/cpp/single/setup/CMakeLists.txt
add_library(setup_single STATIC "setup.cpp")
# set_target_properties(setup_single PROPERTIES OUTPUT_NAME setup)
target_link_libraries(setup_single PUBLIC scatterer source detector)

pybind11_add_module("interface_setup_single" MODULE interface.cpp)
set_target_properties("interface_setup_single" PROPERTIES OUTPUT_NAME setup)
target_link_libraries("interface_setup_single" PUBLIC pybind11::module setup_single)

install(
    TARGETS "setup_single" "interface_setup_single"
    LIBRARY DESTINATION "PyMieSim/single"
    RUNTIME DESTINATION "PyMieSim/single"
    ARCHIVE DESTINATION "PyMieSim/single"
)