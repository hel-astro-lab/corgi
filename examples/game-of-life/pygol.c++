#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "gol.h"
#include "../../toolbox/dataContainer.h"


PYBIND11_MODULE(pygol, m) {


  /// Bind 2D Mesh 
  py::class_<gol::Mesh>(m, "Mesh")
    .def(py::init<int, int>())
    .def_readwrite("Nx",  &gol::Mesh::Nx)
    .def_readwrite("Ny",  &gol::Mesh::Ny)
    .def("__getitem__", [](const gol::Mesh &s, py::tuple indx) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();

        if (i < -s.halo) throw py::index_error();
        if (j < -s.halo) throw py::index_error();

        if (i > s.Nx+s.halo) throw py::index_error();
        if (j > s.Ny+s.halo) throw py::index_error();

        return s(i,j);
      })
    .def("__setitem__", [](gol::Mesh &s, py::tuple indx, int val) 
      {
        int i = indx[0].cast<int>();
        int j = indx[1].cast<int>();

        if (i < -s.halo) throw py::index_error();
        if (j < -s.halo) throw py::index_error();

        if (i > s.Nx+s.halo) throw py::index_error();
        if (j > s.Ny+s.halo) throw py::index_error();

        s(i,j) = val;
      });


  // --------------------------------------------------
  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");

  // cell binding
  py::class_<gol::CellularAutomataCell, 
            corgi::Cell, 
            std::shared_ptr<gol::CellularAutomataCell>
            >(m, "CellularAutomataCell")
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("addData",          &gol::CellularAutomataCell::addData)
    .def("getData",          &gol::CellularAutomataCell::getData)
    .def("cycle",            &gol::CellularAutomataCell::cycle)
    .def("updateBoundaries", &gol::CellularAutomataCell::updateBoundaries);



  // --------------------------------------------------
  // Grid bindings
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<gol::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>());
    // .def("petShop",   &gol::Grid::petShop);

  py::class_<gol::Solver>(m, "Solver")
    .def(py::init<>())
    .def("solve", &gol::Solver::solve);

}


