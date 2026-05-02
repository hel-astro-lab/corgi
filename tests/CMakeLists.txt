# Copyright 2018 - 2026, Joonas Nättilä and the hel-astro-lab contributors
# SPDX-License-Identifier: Apache-2.0

file (GLOB CORGI_TEST_SOURCES *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module (pycorgitest ${CORGI_TEST_SOURCES})
target_link_libraries (pycorgitest PRIVATE corgi corgi_warnings)
