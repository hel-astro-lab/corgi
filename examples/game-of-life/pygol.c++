#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "gol.h"

#include "../../toolbox/dataContainer.h"
#include "../../toolbox/dataContainer.c++"


PYBIND11_MODULE(pygol, m) {


  /// Bind 2D Mesh 
  py::class_<gol::Mesh>(m, "Mesh")
    .def(py::init<size_t, size_t>())
    .def_readwrite("Nx",  &gol::Mesh::Nx)
    .def_readwrite("Ny",  &gol::Mesh::Ny)
    .def("__getitem__", [](const gol::Mesh &s, py::tuple indx) 
      {
        size_t i = indx[0].cast<size_t>();
        size_t j = indx[1].cast<size_t>();

        if (i >= s.Nx) throw py::index_error();
        if (j >= s.Ny) throw py::index_error();

        return s(i,j);
      })
    .def("__setitem__", [](gol::Mesh &s, py::tuple indx, int val) 
      {
        size_t i = indx[0].cast<size_t>();
        size_t j = indx[1].cast<size_t>();

        if (i >= s.Nx) throw py::index_error();
        if (j >= s.Ny) throw py::index_error();

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
    .def("addData", &gol::CellularAutomataCell::addData)
    .def("getData", &gol::CellularAutomataCell::getData);



  // --------------------------------------------------
  // Grid bindings
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<gol::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>());
    // .def("petShop",   &gol::Grid::petShop);


}


