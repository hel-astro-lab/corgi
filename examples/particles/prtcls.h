#pragma once

#include <vector>
#include <mpi4cpp/mpi.h>

#include "corgi/tile.h"
#include "corgi/corgi.h"

#include "container.h"
#include "communication.h"


namespace prtcls {



/// Particle tile
class Tile : public corgi::Tile<2> {

  public:

  using Tile_t = Tile;
  using Tileptr = std::shared_ptr<Tile>;

  Tile() = default;
  ~Tile() override = default;

  /// particle storage
  std::vector<ParticleBlock> containers;

  /// getter
  ParticleBlock& get_container(size_t i) { return containers[i]; };

  /// setter
  void set_container(const ParticleBlock& block) {containers.push_back(block);};

  size_t Nspecies() {return containers.size(); };
    


  //--------------------------------------------------
  // MPI send
  std::vector<mpi4cpp::mpi::request> 
  send_data( mpi4cpp::mpi::communicator& /*comm*/, int dest, int mode, int tag) override;

  /// actual tag=0 send
  std::vector<mpi4cpp::mpi::request> 
  send_particle_data( mpi4cpp::mpi::communicator& /*comm*/, int dest, int tag);

  /// actual tag=1 send
  std::vector<mpi4cpp::mpi::request> 
  send_particle_extra_data( mpi4cpp::mpi::communicator& /*comm*/, int dest, int tag);


  //--------------------------------------------------
  // MPI recv
  std::vector<mpi4cpp::mpi::request> 
  recv_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int mode, int tag) override;

  /// actual tag=0 recv
  std::vector<mpi4cpp::mpi::request> 
  recv_particle_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int tag);

  /// actual tag=1 recv
  std::vector<mpi4cpp::mpi::request> 
  recv_particle_extra_data(mpi4cpp::mpi::communicator& /*comm*/, int orig, int tag);
  //--------------------------------------------------


  /// check all particle containers for particles
  // exceeding limits
  void check_outgoing_particles();

  /// delete particles from each container that are exceeding
  // the boundaries
  void delete_transferred_particles();

  /// get particles flowing into this tile
  void get_incoming_particles(corgi::Grid<2>& grid);

  /// pack particles for MPI message
  void pack_outgoing_particles();

  /// unpack received MPI message particles
  void unpack_incoming_particles();

  /// delete all particles from each container
  void delete_all_particles();

};



/// Particle mover
class Pusher {

  public:
    void solve(Tile& /*tile*/);
};



} // end of namespace prtcls
