file (GLOB CORGI_SOURCES *.c++ *.h)

pybind11_add_module (pycorgi ${CORGI_SOURCES})
target_link_libraries (pycorgi PRIVATE corgi corgi_warnings)

install (TARGETS pycorgi DESTINATION lib)
