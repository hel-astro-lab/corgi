#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "example.h"




PYBIND11_MODULE(example, m) {

  // --------------------------------------------------
  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");

  // py::class_<example::Welsh>(m, "Welsh", corgiCell)
  // py::class_<Derived, Base, std::shared_ptr<Derived>>(...);
  py::class_<example::Welsh, corgi::Cell, std::shared_ptr<example::Welsh>>(m, "Welsh")
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("bark",     &example::Welsh::bark);

  py::class_<example::Pembroke, corgi::Cell, std::shared_ptr<example::Pembroke>>(m, "Pembroke")
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("howl",     &example::Pembroke::howl)
    .def("bark",     &example::Pembroke::bark);


  // --------------------------------------------------
  // Grid bindings
  py::object corgiNode = (py::object) py::module::import("corgi").attr("Node");
  py::class_<example::Grid>(m, "Grid", corgiNode)
    .def(py::init<size_t, size_t>())
    .def("petShop",   &example::Grid::petShop);


}




