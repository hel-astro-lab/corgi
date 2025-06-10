project (pycorgi)

file(GLOB SOURCES *.c++ *.h)

#find_package(pybind11 REQUIRED) #loaded externally

pybind11_add_module(pycorgi ${SOURCES})
target_compile_options(pycorgi PRIVATE ${WARNING_FLAGS})
target_link_libraries(pycorgi PRIVATE MPI::MPI_CXX)


install (TARGETS pycorgi DESTINATION lib)
