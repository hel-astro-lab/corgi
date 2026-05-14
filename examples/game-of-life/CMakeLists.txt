# Copyright 2018 - 2026, Joonas Nättilä and the corgi contributors
# SPDX-License-Identifier: Apache-2.0

file (GLOB CORGI_PYCA_SOURCES *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module (corgi_pyca ${CORGI_PYCA_SOURCES})
target_link_libraries (corgi_pyca PRIVATE corgi corgi_warnings)
