project (cellular_automata)

file (GLOB TEST_SRCS *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module(pyca ${TEST_SRCS})
target_compile_options(pyca PRIVATE ${WARNING_FLAGS})


