#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "gol.h"


PYBIND11_MODULE(example, m) {

  // --------------------------------------------------
  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");

  // cell binding
  py::class_<gol::CellularAutomataCell, 
            corgi::Cell, 
            std::shared_ptr<gol::CellularAutomataCell>
            >(m, "CellularAutomataCell")
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("bark",     &gol::Welsh::bark);


  // --------------------------------------------------
  // Grid bindings
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<gol::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>())
    .def("petShop",   &gol::Grid::petShop);


}


