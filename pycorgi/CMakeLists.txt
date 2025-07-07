file (GLOB CORGI_SOURCES *.c++ *.h)

pybind11_add_module (pycorgi ${CORGI_SOURCES})
target_link_libraries (pycorgi PRIVATE corgi corgi_warnings)

set_target_properties (pycorgi PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/lib"
)
