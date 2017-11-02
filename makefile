

CMP = mpic++
LNK = mpic++

CXXFLAGS+=-Wall -O3 -std=c++14

#pybind in macOS needs to have these additional flags
PYBINDINCLS= `python2 -m pybind11 --includes`
PYBINDFLAGS=-shared -fPIC -undefined dynamic_lookup 



all: pycorgi

default: pycorgi



##################################################
# actual compilation & linking rules

pycorgi/corgi.o: corgi.h common.h toolbox/SparseGrid.h pycorgi/corgi.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o pycorgi/corgi.o -c pycorgi/corgi.c++

pycorgi: pycorgi/corgi.o
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o pycorgi/corgi.so pycorgi/corgi.o


.PHONY: tests
tests: 
	python2 -m unittest discover -s tests/ -v


.PHONY: clean
clean: 
	rm pycorgi/*.o
	rm pycorgi/*.so
