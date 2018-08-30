project (pycorgi)

file(GLOB SOURCES *.c++ *.h)

find_package(pybind11 REQUIRED)

pybind11_add_module(pycorgi ${SOURCES})
target_compile_options(pycorgi PRIVATE ${WARNING_FLAGS})


install (TARGETS pycorgi DESTINATION lib)
