#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

#include "simulation.h"




PYBIND11_MODULE(specialization, m) {

  // --------------------------------------------------
  // Loading cell bindings from corgi library
  py::object corgiCell = (py::object) py::module::import("corgi").attr("Cell");

  py::class_<example::Welsh>(m, "Welsh", corgiCell)
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("bark",     &example::Welsh::bark);

  py::class_<example::Pembroke>(m, "Pembroke", corgiCell)
    .def(py::init<size_t, size_t, int, size_t, size_t>())
    .def("howl",     &example::Pembroke::bark)
    .def("bark",     &example::Pembroke::bark);



  // --------------------------------------------------
  // Grid bindings



}








