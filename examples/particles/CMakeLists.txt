# Copyright 2018 - 2026, Joonas Nättilä and the hel-astro-lab contributors
# SPDX-License-Identifier: Apache-2.0

file (GLOB CORGI_PARTICLES_SOURCES *.c++ *.h *.cpp *.hpp *.c *.cxx)

pybind11_add_module (corgi_pyprtcls ${CORGI_PARTICLES_SOURCES})
target_link_libraries (corgi_pyprtcls PRIVATE corgi corgi_warnings)
