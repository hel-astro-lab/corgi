

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

main.o: ../corgi.h ../common.h ../toolbox/SparseGrid.h main.c++
	${CMP} ${CXXFLAGS} ${PYBINDINCLS} -o main.o -c main.c++

pycorgi: main.o
	${LNK} ${PYBINDFLAGS} ${PYBINDINCLS} ${LDFLAGS} -o corgi.so main.o




.PHONY: clean
clean: 
	rm *.o
	rm *.so
