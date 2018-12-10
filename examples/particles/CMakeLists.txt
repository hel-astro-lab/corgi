project (particles)

file (GLOB TEST_SRCS *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module(pyprtcls ${TEST_SRCS})
target_compile_options(pyprtcls PRIVATE ${WARNING_FLAGS})


