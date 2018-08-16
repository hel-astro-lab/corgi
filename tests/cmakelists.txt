project (pycorgitest)

find_package(pybind11 REQUIRED)

file (GLOB TEST_SRCS *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module(pycorgitest ${TEST_SRCS})
target_compile_options(pycorgitest PRIVATE ${WARNING_FLAGS})
