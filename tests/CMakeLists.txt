file (GLOB CORGI_TEST_SOURCES *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module (pycorgitest ${CORGI_TEST_SOURCES})
target_link_libraries (pycorgitest PRIVATE corgi corgi_warnings)

set_target_properties (pycorgitest PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
)
