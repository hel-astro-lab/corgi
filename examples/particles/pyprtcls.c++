#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;
    
#include "prtcls.h"


PYBIND11_MODULE(pyprtcls, m) {

  // --------------------------------------------------
  py::class_<prtcls::ParticleBlock>(m, "ParticleBlock")
    .def(py::init<size_t, size_t, size_t>())
    .def("reserve",       &prtcls::ParticleBlock::reserve)
    .def("size",          &prtcls::ParticleBlock::size)
    .def("add_particle",  &prtcls::ParticleBlock::add_particle)
    .def("add_particle2", [](prtcls::ParticleBlock& s, 
                            double xx, double yy, double zz,
                            double vx, double vy, double vz, 
                            double wgt)
        {
          s.add_particle({xx,yy,zz}, {vx,vy,vz}, wgt);
        })
    .def("loc",          [](prtcls::ParticleBlock& s, size_t idim) 
        {
          return s.loc(idim); 
        }, py::return_value_policy::reference)
    .def("vel",          [](prtcls::ParticleBlock& s, size_t idim) 
        {
          return s.vel(idim); 
        }, py::return_value_policy::reference)
    .def("wgt",          [](prtcls::ParticleBlock& s) 
        {
          return s.wgt(); 
        }, py::return_value_policy::reference);
    


  // --------------------------------------------------

  // tile binding
  py::class_<prtcls::Tile, 
            corgi::Tile<2>, 
            std::shared_ptr<prtcls::Tile>
            >(m, "Tile")
    .def(py::init<>())
    .def("get_container",       &prtcls::Tile::get_container, 
        py::return_value_policy::reference)
    .def("set_container",       &prtcls::Tile::set_container)
    .def("check_outgoing_particles",     &prtcls::Tile::check_outgoing_particles)
    .def("get_incoming_particles",       &prtcls::Tile::get_incoming_particles)
    .def("delete_transferred_particles", &prtcls::Tile::delete_transferred_particles)
    .def("pack_outgoing_particles",      &prtcls::Tile::pack_outgoing_particles)
    .def("unpack_incoming_particles",    &prtcls::Tile::unpack_incoming_particles)
    .def("delete_all_particles",         &prtcls::Tile::delete_all_particles);


  // --------------------------------------------------
  // Solver

  py::class_<prtcls::Pusher>(m, "Pusher")
    .def(py::init<>())
    .def("solve", &prtcls::Pusher::solve);


}


