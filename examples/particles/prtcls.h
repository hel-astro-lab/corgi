#pragma once

#include <vector>
#include <mpi4cpp/mpi.h>

#include "../../tile.h"
#include "../../corgi.h"

#include "container.h"


namespace prtcls {

namespace mpi = mpi4cpp::mpi;



/// Particle tile
class Tile : public corgi::Tile<2> {

  public:

  typedef Tile Tile_t;
  typedef std::shared_ptr<Tile> Tileptr;

  Tile() {};
  ~Tile() = default;

  /// particle storage
  std::vector<ParticleBlock> containers;

  /// getter
  ParticleBlock& get_container(size_t i) { return containers[i]; };

  /// setter
  void set_container(const ParticleBlock& block) {containers.push_back(block);};


  virtual mpi::request send_data( mpi::communicator&, int orig, int tag) override;
  virtual mpi::request recv_data( mpi::communicator&, int dest, int tag) override;
};



/// Particle mover
class Pusher {

  public:
    void solve(Tile&);
};



} // end of namespace prtcls
